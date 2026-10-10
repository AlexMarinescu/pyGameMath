"""Manual release-candidate validation. There is deliberately no upload operation."""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
REPOSITORY='AlexMarinescu/pyGameMath'
VERSION='1.0.0'
TAG='v1.0.0'
ENVIRONMENT='gem-release-review'
FILES={'gem-1.0.0-py3-none-any.whl','gem-1.0.0.tar.gz'}


def require(condition,message):
    if not condition:raise ValueError(message)


def validate_context(context,commit,version,tag,mode):
    require(context.get('event_name')=='workflow_dispatch','manual dispatch is required')
    require(context.get('repository')==REPOSITORY and not context.get('fork',True),'canonical non-fork repository is required')
    require(context.get('ref')=='refs/heads/master','dispatch must run from master')
    require(re.fullmatch(r'[0-9a-f]{40}',commit or '') is not None,'commit must be a full lowercase SHA')
    require(commit==context.get('sha'),'candidate must match the dispatched master commit')
    require(version==VERSION and tag==TAG,'version/tag must be exactly 1.0.0/v1.0.0')
    require(mode in ('validate','review'),'only validate/review modes are supported; publishing is disabled')
    return {'repository':REPOSITORY,'source_commit':commit,'version':version,'tag':tag,'mode':mode,
            'publication_enabled':False,'ownership_verified':False}


def verify_review_configuration(environment,branch):
    require(branch.get('protected') is True,'master must have server-side branch protection')
    require(environment.get('name')==ENVIRONMENT,'release review environment missing')
    rules=environment.get('protection_rules',[])
    reviewers=[r for r in rules if r.get('type')=='required_reviewers']
    require(len(reviewers)==1 and bool(reviewers[0].get('reviewers')),'environment must require reviewers')
    require(reviewers[0].get('prevent_self_review') is True,'self-review must be prevented')
    require(environment.get('can_admins_bypass') is False,'administrative approval bypass must be disabled')
    policy=environment.get('deployment_branch_policy') or {}
    require(policy.get('protected_branches') is True and policy.get('custom_branch_policies') is False,
            'environment must restrict deployments to protected branches')
    return {'environment':ENVIRONMENT,'required_reviewers':True,'prevent_self_review':True,
            'admin_bypass':False,'protected_branches_only':True}


def github_get(path):
    # GitHub.com is the explicitly supported release repository. Never send a token to an arbitrary host.
    require(os.environ.get('GITHUB_API_URL','https://api.github.com')=='https://api.github.com','unexpected API host')
    token=os.environ.get('GITHUB_TOKEN')
    require(bool(token),'read-only GitHub token is required to inspect review protection')
    request=urllib.request.Request('https://api.github.com/repos/'+REPOSITORY+'/'+path,
        headers={'Authorization':'Bearer '+token,'Accept':'application/vnd.github+json',
                 'X-GitHub-Api-Version':'2022-11-28'})
    with urllib.request.urlopen(request,timeout=30) as response:return json.load(response)


def preflight(context,commit,version,tag,mode,get=github_get):
    report=validate_context(context,commit,version,tag,mode)
    current=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    require(current==commit,'checkout must match candidate commit')
    tree=ast.parse((ROOT/'gem/_version.py').read_text(encoding='utf-8'))
    literal=[ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign)
             and any(isinstance(t,ast.Name) and t.id=='__version__' for t in n.targets)]
    require(literal==[VERSION],'authoritative source version mismatch')
    present=subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag],cwd=ROOT)
    require(present.returncode in (0,1),'cannot determine proposed tag state')
    if present.returncode==0:
        target=subprocess.check_output(['git','rev-parse','--verify','refs/tags/'+tag+'^{commit}'],
                                       cwd=ROOT,text=True).strip()
        require(target==commit,'existing tag points to another commit')
        report['tag_state']='existing tag matches candidate'
    else:report['tag_state']='proposed tag absent; no tag created'
    if mode=='review':
        report['protection']=verify_review_configuration(get('environments/'+ENVIRONMENT),get('branches/master'))
    report['status']='validated; publication remains blocked'
    return report


def verify_artifacts(directory,commit,tag):
    directory=directory.resolve()
    require(set(p.name for p in directory.iterdir())==FILES|{'candidate-manifest.json','SHA256SUMS','INSTALL.md'},'candidate bundle must contain two archives, manifest, checksums and installation instructions')
    require(all(p.is_file() and not p.is_symlink() for p in directory.iterdir()),'candidate files cannot be symlinks/directories')
    manifest=json.loads((directory/'candidate-manifest.json').read_text(encoding='utf-8'))
    require(manifest.get('schema_version')==1 and manifest.get('distribution')=='gem','manifest schema/identity mismatch')
    require(manifest.get('version')==VERSION and manifest.get('tag')==tag==TAG,'artifact version/tag mismatch')
    require(re.fullmatch(r'[0-9a-f]{40}',commit or '') is not None,'full candidate SHA required')
    require(manifest.get('source_commit')==commit,'artifact source commit mismatch')
    require(manifest.get('source_clean') is True,'dirty source cannot become a release candidate')
    tree=subprocess.check_output(['git','rev-parse','HEAD^{tree}'],cwd=ROOT,text=True).strip()
    current=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    require(current==commit and manifest.get('source_tree')==tree,'checkout/tree does not match artifact provenance')
    require(set(manifest.get('files',{}))==FILES,'artifact filename set mismatch')
    for name,digest in manifest['files'].items():
        require(re.fullmatch(r'[0-9a-f]{64}',digest or '') is not None,'invalid SHA-256 digest')
        require(hashlib.sha256((directory/name).read_bytes()).hexdigest()==digest,'artifact digest mismatch: '+name)
    expected=''.join(digest+'  '+name+'\n' for name,digest in sorted(manifest['files'].items()))
    require((directory/'SHA256SUMS').read_text(encoding='utf-8')==expected,'checksum file mismatch')
    require(bool((directory/'INSTALL.md').read_text(encoding='utf-8').strip()),'installation instructions missing')
    spec=importlib.util.spec_from_file_location('distribution_checks',ROOT/'tools/verify_distribution.py')
    checks=importlib.util.module_from_spec(spec);spec.loader.exec_module(checks)
    inspection=checks.inspect_archives(directory/'gem-1.0.0-py3-none-any.whl',directory/'gem-1.0.0.tar.gz')
    return {'status':'artifact integrity verified','manifest':manifest,'runtime_files':inspection['runtime_files'],
            'publication_enabled':False,'ownership_verified':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    gate=sub.add_parser('preflight');gate.add_argument('--output',type=Path,required=True)
    integrity=sub.add_parser('verify-artifacts');integrity.add_argument('--directory',type=Path,required=True)
    integrity.add_argument('--commit',required=True);integrity.add_argument('--tag',required=True)
    integrity.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='preflight':
        event=json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text(encoding='utf-8'))
        context={name:os.environ.get('GITHUB_'+name.upper()) for name in ('event_name','repository','ref','sha')}
        context['fork']=event.get('repository',{}).get('fork',True)
        report=preflight(context,os.environ.get('INPUT_COMMIT'),os.environ.get('INPUT_VERSION'),
                         os.environ.get('INPUT_TAG'),os.environ.get('INPUT_MODE'))
        if os.environ.get('GITHUB_OUTPUT'):
            with open(os.environ['GITHUB_OUTPUT'],'a',encoding='utf-8') as output:
                for name in ('source_commit','version','tag','mode'):output.write(name+'='+report[name]+'\n')
    else:report=verify_artifacts(args.directory,args.commit,args.tag)
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(report['status'])


if __name__=='__main__':main()
