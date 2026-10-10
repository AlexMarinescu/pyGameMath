"""Portable CI orchestration using the Phase 5A checks; never publish."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import sysconfig
import venv
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DOC_FLAGS = ['--examples', '--getting-started', '--api-reference', '--tutorials', '--showcase', '--website']


def environment_python(root, os_name=None):
    return root / ('Scripts/python.exe' if (os_name or os.name) == 'nt' else 'bin/python')


def assert_clean_junit(path):
    suites=ET.parse(path).getroot()
    totals={name:sum(int(s.get(name,0)) for s in suites.iter('testsuite'))
            for name in ('tests','failures','errors','skipped')}
    cases=list(suites.iter('testcase'))
    if totals['tests'] != len(cases) or totals['tests'] == 0 or any(totals[n] for n in ('failures','errors','skipped')) or any(
            child.tag in {'failure','error','skipped'} for case in cases for child in case):
        raise ValueError('full suite must have passing tests and zero failures/errors/skips/xfails')
    return dict(totals, passed=totals['tests'], xfailed=0)


def install_artifacts(output_dir, wheel, sdist, command, read, offline=False, wheelhouse=None):
    installed={}
    for kind,artifact in [('wheel',wheel),('sdist',sdist)]:
        target=output_dir/(kind+'-environment')
        venv.EnvBuilder(with_pip=True,symlinks=os.name!='nt').create(target)
        executable=environment_python(target)
        # Build tools are deliberately separate from gem's runtime dependency.
        setup=[executable,'-m','pip','install','six==1.17.0','setuptools==80.9.0','wheel==0.48.0']
        if offline: setup.extend(['--no-index','--find-links',str(wheelhouse.resolve())])
        command(setup,kind+'-dependencies',cwd=output_dir)
        install=[executable,'-I','-m','pip','install','--no-deps','--no-build-isolation',artifact]
        if offline:install.insert(-1,'--no-index')
        command(install,kind+'-install',cwd=output_dir)
        site=subprocess.check_output([str(executable),'-I','-c',
            'import sysconfig; print(sysconfig.get_path("purelib"))'],text=True,cwd=output_dir).strip()
        command([executable,'-I','-X','utf8',ROOT/'tools/verify_distribution.py','--installed-root',site,
                 '--output',output_dir/(kind+'-smoke.json')],kind+'-smoke',cwd=output_dir)
        command([executable,'-I','-X','utf8',ROOT/'tools/check_architecture_docs.py',*DOC_FLAGS,'--package-root',site,
                 '--output',output_dir/(kind+'-docs.json')],kind+'-docs',cwd=output_dir)
        docs=read(kind+'-docs.json')
        if docs['source_declarations_checked']!=268 or docs['examples_executed']!=43:
            raise ValueError('installed API/example coverage incomplete')
        installed[kind]={'smoke':read(kind+'-smoke.json'),'documentation':docs}
    return installed


def run(output_dir, artifacts_dir, offline=False, wheelhouse=None):
    if sys.flags.optimize: raise ValueError('CI assertions must not run with Python optimization enabled')
    if offline and wheelhouse is None: raise ValueError('offline verification requires an explicit wheelhouse')
    output_dir=output_dir.resolve();artifacts_dir=artifacts_dir.resolve()
    output_dir.mkdir(parents=True,exist_ok=True)
    if artifacts_dir.exists() and any(artifacts_dir.iterdir()):
        raise ValueError('artifact directory must be empty for a clean build')
    artifacts_dir.mkdir(parents=True,exist_ok=True)
    initially_clean=not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    report={'python':platform.python_version(),'platform':platform.platform(),'commands':[],
            'status':'running','working_directory':str(ROOT)}
    env=os.environ.copy();env.pop('PYTHONPATH',None);env['PIP_DISABLE_PIP_VERSION_CHECK']='1';env['PYTHONUTF8']='1'
    def command(args, name, cwd=ROOT):
        report['commands'].append({'argv':list(map(str,args)),'cwd':str(cwd),'log':name+'.log'})
        with (output_dir/(name+'.log')).open('w',encoding='utf-8') as log:
            subprocess.run(list(map(str,args)),cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    def read(name):return json.loads((output_dir/name).read_text(encoding='utf-8'))
    try:
        command([sys.executable,ROOT/'tools/trace_showcase.py','--output',
                 output_dir/'showcase-trace.json'],'showcase-trace')
        pytest_error=None
        try:
            command([sys.executable,'-m','pytest','-q','-o','junit_family=legacy',
                     '--junitxml',str(output_dir/'pytest.xml')],'pytest')
            report['tests']=assert_clean_junit(output_dir/'pytest.xml')
        except (subprocess.CalledProcessError,ValueError) as error:
            # Continue independent archive/install diagnostics, but never approve a failed suite.
            pytest_error=error
            suites=ET.parse(output_dir/'pytest.xml').getroot()
            totals={name:sum(int(s.get(name,0)) for s in suites.iter('testsuite'))
                    for name in ('tests','failures','errors','skipped')}
            report['tests']=dict(totals,passed=totals['tests']-sum(totals[n] for n in ('failures','errors','skipped')))
        command([sys.executable,ROOT/'tools/check_architecture_docs.py',*DOC_FLAGS,
                 '--output',output_dir/'source-docs.json'],'source-docs')
        report['documentation']=read('source-docs.json')
        if report['documentation']['source_declarations_checked']!=268 or report['documentation']['examples_executed']!=43:
            raise ValueError('API/example coverage must remain 268/43')
        build=[sys.executable,'-m','build','--outdir',artifacts_dir]
        if offline:build.append('--no-isolation')
        command(build,'build')
        wheel=artifacts_dir/'gem-1.0.0-py3-none-any.whl';sdist=artifacts_dir/'gem-1.0.0.tar.gz'
        if set(p.name for p in artifacts_dir.iterdir())!={wheel.name,sdist.name}:
            raise ValueError('build must produce exactly the expected wheel and sdist')
        command([sys.executable,'-m','twine','check','--strict',wheel,sdist],'twine')
        command([sys.executable,ROOT/'tools/verify_distribution.py','--wheel',wheel,'--sdist',sdist,
                 '--output',output_dir/'artifacts.json'],'artifacts')
        report['artifacts']=read('artifacts.json')
        report['installed']=install_artifacts(output_dir,wheel,sdist,command,read,offline,wheelhouse)
        if pytest_error is not None:
            raise ValueError('suite failed; archive/install diagnostics do not approve the candidate') from pytest_error
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
        tree=subprocess.check_output(['git','rev-parse','HEAD^{tree}'],cwd=ROOT,text=True).strip()
        manifest={'schema_version':1,'distribution':'gem','version':'1.0.0','tag':'v1.0.0',
                  'source_commit':commit,'source_tree':tree,
                  'source_clean':initially_clean and not subprocess.check_output(
                      ['git','status','--porcelain'],cwd=ROOT,text=True).strip(),
                  'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (wheel,sdist)}}
        (artifacts_dir/'SHA256SUMS').write_text(''.join(
            digest+'  '+name+'\n' for name,digest in sorted(manifest['files'].items())),encoding='utf-8')
        (artifacts_dir/'INSTALL.md').write_text(
            '# gem 1.0.0 candidate installation\n\n'
            'Verify SHA256SUMS against the approved release checksums before installation.\n'
            'Use CPython 3.10–3.14; six is the sole runtime dependency.\n\n'
            'python -m pip install ./gem-1.0.0-py3-none-any.whl\n\n'
            'Alternatively: python -m pip install ./gem-1.0.0.tar.gz\n\n'
            'These are validated candidate assets, not evidence of publication.\n'
            'GitHub-first release publication requires explicit owner approval.\n'
            'PyPI publishing authority remains unverified. Preserve the approved archives\n'
            'for later PyPI publication; never silently rebuild official artifacts.\n',encoding='utf-8')
        (artifacts_dir/'candidate-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
        report['candidate_manifest']=manifest;report['status']='passed'
    except Exception as error:
        report['status']='failed';report['error_type']=type(error).__name__
        raise
    finally:
        (output_dir/'ci-results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report


def verify_candidate(output_dir, directory, offline=False, wheelhouse=None):
    """Install the same canonical archives on each platform without rebuilding."""
    from release_gate import verify_artifacts
    if sys.flags.optimize: raise ValueError('CI assertions must not run with Python optimization enabled')
    if offline and wheelhouse is None: raise ValueError('offline verification requires an explicit wheelhouse')
    output_dir=output_dir.resolve(); directory=directory.resolve()
    output_dir.mkdir(parents=True,exist_ok=True)
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    report={'python':platform.python_version(),'platform':platform.platform(),'commands':[],
            'status':'running','source_commit':commit}
    env=os.environ.copy();env.pop('PYTHONPATH',None);env['PYTHONUTF8']='1'
    env['PIP_DISABLE_PIP_VERSION_CHECK']='1'
    def command(args, name, cwd=ROOT):
        report['commands'].append({'argv':list(map(str,args)),'cwd':str(cwd),'log':name+'.log'})
        with (output_dir/(name+'.log')).open('w',encoding='utf-8') as log:
            subprocess.run(list(map(str,args)),cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    def read(name):return json.loads((output_dir/name).read_text(encoding='utf-8'))
    try:
        report['integrity']=verify_artifacts(directory,commit,'v1.0.0')
        report['installed']=install_artifacts(output_dir,directory/'gem-1.0.0-py3-none-any.whl',
            directory/'gem-1.0.0.tar.gz',command,read,offline,wheelhouse)
        report['status']='passed'
    except Exception as error:
        report['status']='failed';report['error_type']=type(error).__name__
        raise
    finally:
        (output_dir/'candidate-install-results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--artifacts-dir',type=Path)
    parser.add_argument('--candidate-directory',type=Path,help='verify/install an existing canonical bundle; never rebuild')
    parser.add_argument('--offline',action='store_true',help='use preinstalled pinned dependencies; no index/build isolation')
    parser.add_argument('--wheelhouse',type=Path,help='local wheels for offline installed dependencies')
    args=parser.parse_args()
    if bool(args.artifacts_dir) == bool(args.candidate_directory):
        parser.error('choose exactly one of --artifacts-dir or --candidate-directory')
    report=(verify_candidate(args.output_dir,args.candidate_directory,args.offline,args.wheelhouse)
            if args.candidate_directory else run(args.output_dir,args.artifacts_dir,args.offline,args.wheelhouse))
    print(json.dumps({'status':report['status'],'python':report['python'],'tests':report.get('tests')},indent=2))


if __name__=='__main__':main()
