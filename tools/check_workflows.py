"""Validate repository-specific CI safety invariants (PyYAML is a CI-only tool)."""
import argparse
import copy
import json
from pathlib import Path
import re

import yaml

ROOT=Path(__file__).resolve().parents[1]


class WorkflowLoader(yaml.SafeLoader):
    pass


# YAML 1.2 keeps the workflow key `on` as text, unlike PyYAML's default YAML 1.1.
WorkflowLoader.yaml_implicit_resolvers=copy.deepcopy(yaml.SafeLoader.yaml_implicit_resolvers)
for key,rules in WorkflowLoader.yaml_implicit_resolvers.items():
    WorkflowLoader.yaml_implicit_resolvers[key]=[(tag,regex) for tag,regex in rules
                                                if tag!='tag:yaml.org,2002:bool']
WorkflowLoader.add_implicit_resolver('tag:yaml.org,2002:bool',re.compile(r'^(?:true|false|True|False|TRUE|FALSE)$'),list('tTfF'))


def load_workflows(root=ROOT):
    return {p.name:yaml.load(p.read_text(encoding='utf-8'),Loader=WorkflowLoader)
            for p in (root/'.github/workflows').glob('*.yml')}


def require(condition,message):
    if not condition:raise ValueError(message)


def validate(workflows):
    expected={'package-validation.yml','documentation-artifact.yml','publish-docs.yml','release-validation.yml','github-release.yml'}
    require(set(workflows)==expected,'workflow inventory differs from reviewed CI baseline')
    pins=[]
    for name,workflow in workflows.items():
        require(workflow.get('permissions')=={'contents':'read'},'default permissions must be read-only: '+name)
        require(not {'pull_request_target','workflow_run','release'} & set(workflow['on']),'unsafe/unreviewed event: '+name)
        for job_name,job in workflow['jobs'].items():
            require(not job.get('continue-on-error'),'failure gate bypass: '+name)
            permissions=job.get('permissions',workflow['permissions'])
            expected_permissions={'contents':'read'}
            if name=='publish-docs.yml' and job_name=='deploy':
                expected_permissions.update({'pages':'write','id-token':'write'})
            elif name=='release-validation.yml' and job_name=='preflight':
                expected_permissions['actions']='read'
            elif name=='github-release.yml':
                expected_permissions['actions']='read'
                if job_name=='publish':expected_permissions['contents']='write'
            require(permissions==expected_permissions,'unexpected privileged permissions: '+name+'/'+job_name)
            require(job.get('secrets')!='inherit','secrets inheritance is forbidden')
            if 'uses' in job:
                require(job['uses'] in {'./.github/workflows/package-validation.yml',
                    './.github/workflows/documentation-artifact.yml'},'unreviewed reusable workflow')
            for step in job.get('steps',[]):
                require(not step.get('continue-on-error'),'step failure bypass')
                if 'uses' in step:
                    require(re.fullmatch(r'[\w-]+/[\w-]+@[0-9a-f]{40}',step['uses']) is not None,'actions must use full commit pins')
                    require('pypi-publish' not in step['uses'],'package publication must remain disabled')
                    pins.append(step['uses'])
                    if step['uses'].startswith('actions/checkout@'):
                        require(step.get('with',{}).get('persist-credentials') is False,'checkout credentials must not persist')
                command=step.get('run','')
                require(not re.search(r'twine\s+upload|gh\s+release\s+create|git\s+(?:push|tag)\b|prepare_wiki',command),
                        'publication/tag/wiki mutation command in active workflow')
                require('secrets.' not in str(step),'no publishing secret wiring is permitted')
    package=workflows['package-validation.yml']
    require(set(package['on'])=={'pull_request','push','workflow_dispatch','workflow_call'},'package events changed')
    require(package['on']['push'].get('branches')==['master'],'package push must target master')
    matrix=package['jobs']['package']['strategy']['matrix']
    require(matrix=={'os':['ubuntu-latest','windows-latest','macos-latest'],
                     'python':['3.10','3.11','3.12','3.13','3.14']},'complete Python/platform matrix is required')
    require(package['jobs']['package']['strategy'].get('fail-fast') is False,'all platforms must report')
    require(any('tools/ci_validate.py' in s.get('run','') for s in package['jobs']['package']['steps']),
            'shared packaging validation missing')
    canonical=package['jobs'].get('candidate-install',{})
    require(canonical.get('needs')=='package','canonical installs must await whole matrix')
    require(canonical.get('strategy',{}).get('matrix')==matrix and
            canonical['strategy'].get('fail-fast') is False,'canonical installation matrix incomplete')
    downloads=[s for s in canonical.get('steps',[]) if s.get('uses','').startswith('actions/download-artifact@')]
    require(len(downloads)==1 and downloads[0].get('with',{}).get('digest-mismatch')=='error' and
            downloads[0]['with'].get('name')=='gem-candidate-${{ github.run_id }}',
            'canonical installs require verified same-run artifacts')
    require(any('--candidate-directory' in s.get('run','') for s in canonical.get('steps',[])),
            'canonical archive verification missing')
    docs=workflows['documentation-artifact.yml'];steps=docs['jobs']['build']['steps']
    commands=[s.get('run','') for s in steps]
    require(any('check_architecture_docs.py' in c and '--api-reference' in c and '--tutorials' in c for c in commands),
            'API/tutorial checks missing')
    require(any('mkdocs build --strict' in c for c in commands) and any('check_site.py' in c for c in commands),
            'strict emitted-site checks missing')
    upload=next(i for i,s in enumerate(steps) if s.get('uses','').startswith('actions/upload-pages-artifact@'))
    require(all(i<upload for i,s in enumerate(steps) if any(x in s.get('run','') for x in ('check_architecture_docs','mkdocs build','check_site'))),
            'Pages upload precedes checks')
    require('github.ref' in steps[upload].get('if','') and 'github.repository' in steps[upload].get('if',''),
            'Pages artifact must be master/canonical-only')
    publish=workflows['publish-docs.yml']
    require(set(publish['on'])=={'push','workflow_dispatch'} and publish['on']['push']['branches']==['master'],
            'Pages publication events changed')
    require(publish['jobs']['deploy']['needs']=='build' and
            publish['jobs']['build']['uses']=='./.github/workflows/documentation-artifact.yml',
            'Pages deployment must reuse completed validation')
    release=workflows['release-validation.yml']
    require(set(release['on'])=={'workflow_dispatch'},'release workflow must be dispatch-only')
    inputs=release['on']['workflow_dispatch']['inputs']
    require(inputs['mode']['options']==['validate','review'] and inputs['mode']['default']=='validate',
            'publication mode is forbidden')
    require(inputs['version']['default']=='1.0.0' and inputs['tag']['default']=='v1.0.0','exact version/tag defaults required')
    guard="github.repository == 'AlexMarinescu/pyGameMath' && github.ref == 'refs/heads/master' && github.event_name == 'workflow_dispatch'"
    jobs=release['jobs'];require(jobs['preflight']['if']==guard,'release source/event guard changed')
    require(jobs['validation']['needs']=='preflight' and jobs['validation']['uses']=='./.github/workflows/package-validation.yml',
            'release must reuse packaging checks after preflight')
    require(set(jobs['integrity']['needs'])=={'preflight','validation'},'artifact checks must await whole matrix')
    require(set(jobs['review']['needs'])=={'preflight','validation','integrity'} and
            jobs['review']['environment']=='gem-release-review' and jobs['review']['if']=="needs.preflight.outputs.mode == 'review'",
            'human review gate bypass')
    github_release=workflows['github-release.yml']
    require(set(github_release['on'])=={'workflow_dispatch'},'GitHub release must be manual-only')
    gjobs=github_release['jobs']
    require(gjobs['validate']['if']==guard,'GitHub release source/event guard changed')
    require(gjobs['publish']['needs']=='validate' and gjobs['publish']['environment']=='gem-release-review',
            'GitHub release protected approval bypass')
    require(gjobs['publish']['if']=="inputs.mode == 'publish' && vars.GEM_RELEASE_ENABLED == 'true' && inputs.owner_authorization == 'publish v1.0.0'",
            'GitHub publication requires explicit enablement and authorization')
    require(github_release['on']['workflow_dispatch']['inputs']['mode']['default']=='validate',
            'GitHub release must default to dry validation')
    return {'workflows':sorted(workflows),'action_pins':sorted(set(pins)),
            'matrix_jobs':15,'canonical_install_jobs':15,'package_publication_enabled':False,'wiki_mutation_enabled':False,
            'pages_write_scope':'master-only deployment job','manual_release_modes':['validate','review']}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=validate(load_workflows());text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text,encoding='utf-8')
    print(text)


if __name__=='__main__':main()
