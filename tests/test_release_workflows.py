"""Independent negative checks for CI and release-candidate boundaries."""
import copy
import importlib.util
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'tools' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gate = load('release_gate')
ci = load('ci_validate')
archives = load('verify_distribution')
workflows = load('check_workflows')
SHA = 'a' * 40
CONTEXT = dict(event_name='workflow_dispatch', repository=gate.REPOSITORY,
               fork=False, ref='refs/heads/master', sha=SHA)


def test_manual_validation_never_enables_publication():
    for mode in ('validate', 'review'):
        result = gate.validate_context(CONTEXT, SHA, '1.0.0', 'v1.0.0', mode)
        assert result['publication_enabled'] is False
        assert result['ownership_verified'] is False


@pytest.mark.parametrize('field,value', [
    ('event_name', 'push'), ('event_name', 'pull_request'),
    ('repository', 'someone/pyGameMath'), ('fork', True),
    ('ref', 'refs/heads/topic'), ('sha', 'b' * 40)])
def test_reject_unreviewed_context(field, value):
    context = dict(CONTEXT, **{field: value})
    with pytest.raises(ValueError):
        gate.validate_context(context, SHA, '1.0.0', 'v1.0.0', 'validate')


@pytest.mark.parametrize('commit,version,tag,mode', [
    ('a' * 7, '1.0.0', 'v1.0.0', 'validate'),
    (SHA.upper(), '1.0.0', 'v1.0.0', 'validate'),
    (SHA, '1.0.1', 'v1.0.0', 'validate'),
    (SHA, '1.0.0', 'v1.0.1', 'validate'),
    (SHA, '1.0.0', 'v1.0.0', 'publish')])
def test_reject_unapproved_candidate(commit, version, tag, mode):
    with pytest.raises(ValueError):
        gate.validate_context(CONTEXT, commit, version, tag, mode)


def protected_environment():
    return dict(name=gate.ENVIRONMENT, can_admins_bypass=False,
                protection_rules=[dict(type='required_reviewers', reviewers=[{'id': 1}],
                                       prevent_self_review=True)],
                deployment_branch_policy=dict(protected_branches=True,
                                              custom_branch_policies=False))


@pytest.mark.parametrize('defect', ['reviewers', 'self_review', 'bypass', 'branches', 'protection'])
def test_review_requires_server_side_protection(defect):
    environment = protected_environment()
    branch = {'protected': True}
    assert gate.verify_review_configuration(environment, branch)['required_reviewers']
    if defect == 'reviewers': environment['protection_rules'][0]['reviewers'] = []
    if defect == 'self_review': environment['protection_rules'][0]['prevent_self_review'] = False
    if defect == 'bypass': environment['can_admins_bypass'] = True
    if defect == 'branches': environment['deployment_branch_policy']['protected_branches'] = False
    if defect == 'protection': branch['protected'] = False
    with pytest.raises(ValueError):
        gate.verify_review_configuration(environment, branch)


@pytest.mark.parametrize('returncode', [0, 1, 2])
def test_tag_detection_fails_closed(monkeypatch, returncode):
    monkeypatch.setattr(gate.subprocess, 'check_output', lambda *a, **k: SHA + '\n')
    monkeypatch.setattr(gate.subprocess, 'run', lambda *a, **k: subprocess.CompletedProcess(a, returncode))
    if returncode == 2:
        with pytest.raises(ValueError): gate.preflight(CONTEXT, SHA, '1.0.0', 'v1.0.0', 'validate')
    else:
        result = gate.preflight(CONTEXT, SHA, '1.0.0', 'v1.0.0', 'validate')
        assert ('absent' in result['tag_state']) == (returncode == 1)


def test_cross_platform_environment_paths(tmp_path):
    assert ci.environment_python(tmp_path, 'nt') == tmp_path / 'Scripts/python.exe'
    assert ci.environment_python(tmp_path, 'posix') == tmp_path / 'bin/python'


@pytest.mark.parametrize('child,tests', [('', 1), ('<failure/>', 1), ('<error/>', 1),
                                        ('<skipped type="pytest.xfail"/>', 1), ('', 0), ('', 2)])
def test_junit_gate_checks_actual_cases(tmp_path, child, tests):
    path = tmp_path / 'results.xml'
    path.write_text('<testsuites><testsuite tests="%d"><testcase>%s</testcase></testsuite></testsuites>' % (tests, child))
    if tests == 1 and not child:
        assert ci.assert_clean_junit(path)['passed'] == 1
    else:
        with pytest.raises(ValueError): ci.assert_clean_junit(path)


@pytest.mark.parametrize('name', ['../secret', '/absolute', 'gem/../secret', 'gem//bad',
                                 'gem\\bad', '.git/config', 'gem/__pycache__/a.pyc',
                                 'dist/a.whl', '.env', 'secret.pem', 'audit/cache.tmp'])
def test_unsafe_archive_members(name):
    with pytest.raises(ValueError): archives.check_member(name)


def test_credential_guard_does_not_echo_credentials():
    token = b'pypi-' + b'x' * 60
    with pytest.raises(ValueError) as error: archives.check_member('example.txt', token)
    assert token.decode() not in str(error.value)
    archives.check_member('tools/site/check_site.py', b'ordinary source')


def test_reviewed_workflow_inventory():
    result = workflows.validate(workflows.load_workflows())
    assert result['matrix_jobs'] == 15
    assert result['package_publication_enabled'] is False


@pytest.mark.parametrize('defect', ['permissions', 'event', 'pin', 'mode', 'review', 'matrix'])
def test_workflow_policy_rejects_bypasses(defect):
    data = copy.deepcopy(workflows.load_workflows())
    release = data['release-validation.yml']
    if defect == 'permissions': release['permissions']['contents'] = 'write'
    if defect == 'event': release['on']['push'] = None
    if defect == 'pin': release['jobs']['preflight']['steps'][0]['uses'] = 'actions/checkout@main'
    if defect == 'mode': release['on']['workflow_dispatch']['inputs']['mode']['options'].append('publish')
    if defect == 'review': release['jobs']['review']['environment'] = 'unprotected'
    if defect == 'matrix': data['package-validation.yml']['jobs']['package']['strategy']['matrix']['python'].pop()
    with pytest.raises(ValueError): workflows.validate(data)


release = load('github_release')


@pytest.mark.parametrize('enabled,approved', [('false', 'publish v1.0.0'),
                                             ('true', ''), (None, None)])
def test_github_publication_requires_explicit_enablement(enabled, approved):
    with pytest.raises(ValueError):
        release.validate_publication(CONTEXT, enabled, approved, 'publish')
    release.validate_publication(CONTEXT, enabled, approved, 'validate')


def test_github_publication_requires_canonical_context():
    release.validate_publication(CONTEXT, 'true', 'publish v1.0.0', 'publish')
    with pytest.raises(ValueError):
        release.validate_publication(dict(CONTEXT, fork=True), 'true', 'publish v1.0.0', 'publish')


@pytest.mark.parametrize('field,value', [('head_sha', 'b' * 40), ('event', 'push'),
                                        ('head_branch', 'topic'), ('conclusion', 'failure'),
                                        ('status', 'in_progress'), ('path', 'other.yml')])
def test_github_release_requires_exact_successful_candidate_run(field, value):
    run = dict(id=123, repository={'full_name': gate.REPOSITORY},
               path='.github/workflows/release-validation.yml', event='workflow_dispatch',
               head_branch='master', head_sha=SHA, status='completed', conclusion='success')
    release.validate_run(run, SHA, '123')
    run[field] = value
    with pytest.raises(ValueError): release.validate_run(run, SHA, '123')


def test_github_publication_guard_cannot_be_removed():
    data = workflows.load_workflows()
    data['github-release.yml']['jobs']['publish']['if'] = "inputs.mode == 'publish'"
    with pytest.raises(ValueError): workflows.validate(data)


@pytest.mark.parametrize('mode', ['validate', 'publish'])
def test_release_execution_uses_exact_assets_and_no_dry_run_mutation(tmp_path, monkeypatch, mode):
    import json
    import urllib.error
    event = tmp_path / 'event.json'
    event.write_text(json.dumps({'repository': {'fork': False}}))
    for key, value in {'GITHUB_EVENT_PATH': str(event), 'GITHUB_EVENT_NAME': 'workflow_dispatch',
                       'GITHUB_REPOSITORY': gate.REPOSITORY, 'GITHUB_REF': 'refs/heads/master',
                       'GITHUB_SHA': SHA, 'RELEASE_ENABLED': 'true',
                       'OWNER_AUTHORIZATION': 'publish v1.0.0'}.items():
        monkeypatch.setenv(key, value)
    run = dict(id=123, repository={'full_name': gate.REPOSITORY},
               path='.github/workflows/release-validation.yml', event='workflow_dispatch',
               head_branch='master', head_sha=SHA, status='completed', conclusion='success')
    def get(path):
        if path == 'actions/runs/123': return run
        if path.startswith('environments/'): return protected_environment()
        if path == 'branches/master': return {'protected': True}
        if path == 'git/ref/tags/v1.0.0' and commands:
            return {'object': {'type': 'commit', 'sha': SHA}}
        raise urllib.error.HTTPError('https://api.github.com', 404, 'absent', {}, None)
    monkeypatch.setattr(release.gate, 'github_get', get)
    preflights = []
    monkeypatch.setattr(release.gate, 'preflight', lambda *a: preflights.append(a))
    monkeypatch.setattr(release.gate, 'verify_artifacts', lambda *a: {'publication_enabled': False})
    commands = []
    monkeypatch.setattr(release.subprocess, 'run', lambda argv, **kw: commands.append(argv))
    monkeypatch.setattr('sys.argv', ['github_release.py', '--directory', str(tmp_path),
                                  '--commit', SHA, '--run-id', '123', '--mode', mode,
                                  '--output', str(tmp_path / 'report.json')])
    release.main()
    assert len(preflights) == 1
    if mode == 'validate':
        assert commands == []
    else:
        assert len(commands) == 2
        assert commands[0][-2:] == ['-f', 'sha=' + SHA]
        command = commands[-1]
        assert command[:4] == ['gh', 'release', 'create', 'v1.0.0']
        assert command[command.index('--target') + 1] == SHA
        assert set(command[-4:]) == {str(tmp_path / n) for n in gate.FILES | {'SHA256SUMS', 'INSTALL.md'}}
        assert '--clobber' not in command and '--verify-tag' in command


def test_ci_scope_cannot_exempt_mathematics(tmp_path, monkeypatch):
    checker = load('check_architecture_docs')
    import json
    (tmp_path / 'tools').mkdir()
    (tmp_path / 'tools/phase5b-ci-scope.json').write_text(json.dumps(
        {'base': checker.BASE, 'files': {'gem/vector.py': '0' * 64}}))
    monkeypatch.setattr(checker, 'ROOT', tmp_path)
    monkeypatch.setattr(checker.subprocess, 'check_output', lambda *a, **k: '')
    with pytest.raises(ValueError, match='mathematics cannot be exempted'):
        checker.check_scope()


@pytest.mark.parametrize('kind', ['commit', 'tag'])
def test_remote_release_tag_must_resolve_to_approved_commit(monkeypatch, kind):
    def get(path):
        return {'object': {'type': 'commit' if path.startswith('git/tags/') else kind,
                           'sha': 'b' * 40}}
    monkeypatch.setattr(release.gate, 'github_get', get)
    with pytest.raises(ValueError, match='another commit'):
        release.remote_tag_matches(SHA)


@pytest.mark.parametrize('case,reason', [
    ('extra_file', 'candidate bundle'), ('dirty_source', 'dirty source'),
    ('wrong_commit', 'source commit'), ('wrong_tree', 'checkout/tree'),
    ('wrong_version', 'version/tag'), ('digest', 'digest mismatch'),
    ('checksums', 'checksum file'), ('empty_instructions', 'instructions missing')])
def test_candidate_integrity_rejects_tampering(tmp_path, monkeypatch, case, reason):
    import hashlib
    import json
    tree = 'c' * 40
    files = {name: hashlib.sha256(name.encode()).hexdigest() for name in gate.FILES}
    manifest = dict(schema_version=1, distribution='gem', version='1.0.0', tag='v1.0.0',
                    source_commit=SHA, source_tree=tree, source_clean=True, files=files)
    for name in files: (tmp_path / name).write_bytes(name.encode())
    (tmp_path / 'SHA256SUMS').write_text(''.join(d + '  ' + n + '\n' for n, d in sorted(files.items())))
    (tmp_path / 'INSTALL.md').write_text('Installation instructions')
    if case == 'extra_file': (tmp_path / 'extra.txt').write_text('unexpected')
    if case == 'dirty_source': manifest['source_clean'] = False
    if case == 'wrong_commit': manifest['source_commit'] = 'b' * 40
    if case == 'wrong_tree': manifest['source_tree'] = 'b' * 40
    if case == 'wrong_version': manifest['version'] = '1.0.1'
    if case == 'digest': files[next(iter(files))] = '0' * 64
    if case == 'checksums': (tmp_path / 'SHA256SUMS').write_text('wrong')
    if case == 'empty_instructions': (tmp_path / 'INSTALL.md').write_text('')
    (tmp_path / 'candidate-manifest.json').write_text(json.dumps(manifest))
    monkeypatch.setattr(gate.subprocess, 'check_output',
                        lambda argv, **kwargs: (tree if argv[-1] == 'HEAD^{tree}' else SHA) + '\n')
    with pytest.raises(ValueError, match=reason):
        gate.verify_artifacts(tmp_path, SHA, 'v1.0.0')


def test_candidate_scope_cannot_exempt_runtime(tmp_path, monkeypatch):
    import json
    checker = load('check_architecture_docs')
    monkeypatch.setattr(checker, 'ROOT', tmp_path)
    monkeypatch.setattr(checker.subprocess, 'check_output', lambda *args, **kwargs: '')
    (tmp_path/'tools').mkdir()
    (tmp_path/'tools/phase5c-reference-scope.json').write_text(json.dumps({
        'base': checker.BASE, 'files': {'gem/quaternion.py': 'unauthorized'}}))
    with pytest.raises(ValueError, match='mathematics cannot be exempted'):
        checker.check_scope()
