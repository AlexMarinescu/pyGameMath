"""Validate an approved GitHub candidate; mutation requires an explicit gated invocation."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import urllib.error

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('release_gate', ROOT / 'tools/release_gate.py')
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def validate_run(run, commit, run_id):
    gate.require(str(run.get('id')) == str(run_id), 'validation run ID mismatch')
    gate.require(run.get('repository', {}).get('full_name') == gate.REPOSITORY,
                 'candidate run must belong to canonical repository')
    gate.require(run.get('path') == '.github/workflows/release-validation.yml', 'wrong validation workflow')
    gate.require(run.get('event') == 'workflow_dispatch' and run.get('head_branch') == 'master',
                 'candidate must come from a manual master validation')
    gate.require(run.get('head_sha') == commit, 'candidate run source mismatch')
    gate.require(run.get('status') == 'completed' and run.get('conclusion') == 'success',
                 'candidate validation must be completed successfully')


def validate_publication(context, enabled, approved, mode):
    gate.validate_context(context, context.get('sha'), gate.VERSION, gate.TAG, 'validate')
    gate.require(mode in ('validate', 'publish'), 'unknown GitHub release mode')
    if mode == 'publish':
        gate.require(enabled == 'true', 'GitHub publication switch is disabled')
        gate.require(approved == 'publish v1.0.0', 'exact owner authorization phrase required')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--commit', required=True)
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--mode', choices=['validate', 'publish'], default='validate')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    event = json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text(encoding='utf-8'))
    context = {name: os.environ.get('GITHUB_' + name.upper())
               for name in ('event_name', 'repository', 'ref', 'sha')}
    context['fork'] = event.get('repository', {}).get('fork', True)
    validate_publication(context, os.environ.get('RELEASE_ENABLED'),
                         os.environ.get('OWNER_AUTHORIZATION'), args.mode)
    gate.require(args.commit == context['sha'], 'release must match final dispatched master commit')
    gate.require(args.run_id.isdecimal(), 'numeric validation run ID required')
    run = gate.github_get('actions/runs/' + args.run_id)
    validate_run(run, args.commit, args.run_id)
    protection = gate.verify_review_configuration(
        gate.github_get('environments/' + gate.ENVIRONMENT), gate.github_get('branches/master'))
    gate.preflight(context, args.commit, gate.VERSION, gate.TAG, 'review')
    result = gate.verify_artifacts(args.directory, args.commit, gate.TAG)
    result.update(mode=args.mode, validation_run_id=args.run_id, protection=protection,
                  github_publication_enabled=args.mode == 'publish', pypi_publication_enabled=False)
    # A release is immutable by policy: never overwrite assets or reuse a published release.
    try:
        gate.github_get('releases/tags/' + gate.TAG)
    except urllib.error.HTTPError as error:
        if error.code != 404: raise
    else:
        raise ValueError('release already exists; refuse to overwrite official artifacts')
    if args.mode == 'publish':
        # --target binds tag creation to the reviewed SHA; existing tags were checked by preflight.
        subprocess.run(['gh', 'release', 'create', gate.TAG, '--repo', gate.REPOSITORY,
                        '--target', args.commit, '--title', 'gem 1.0.0',
                        '--notes-file', str(ROOT / 'docs/development/release-notes-1.0.0.md'),
                        *[str(args.directory / name) for name in sorted(gate.FILES | {'SHA256SUMS', 'INSTALL.md'})]],
                       check=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
