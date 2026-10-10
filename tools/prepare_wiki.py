"""Guarded local preparation of reviewed Wiki pages; never commits or pushes."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT/'docs/wiki-migration'


def prepare(wiki_dir, apply=False):
    wiki_dir = Path(wiki_dir).resolve()
    manifest = json.loads((PACKAGE/'manifest.json').read_text())
    tip = subprocess.check_output(['git', '-C', str(wiki_dir), 'rev-parse', 'HEAD'], text=True).strip()
    if tip != manifest['observed_tip']:
        raise ValueError('Wiki changed since inspection; refresh and review the migration')
    if subprocess.check_output(['git', '-C', str(wiki_dir), 'status', '--porcelain'], text=True).strip():
        raise ValueError('Wiki worktree must be clean')
    observed = manifest.get('observed_files')
    if observed is None:
        if (wiki_dir/'Home.md').read_bytes() != (PACKAGE/'pages/Historical-Home.md').read_bytes():
            raise ValueError('historical Home differs from preserved source')
    else:
        for name, digest in observed.items():
            if Path(name).name != name or hashlib.sha256((wiki_dir/name).read_bytes()).hexdigest() != digest:
                raise ValueError('observed Wiki file differs: '+name)
    # Validate every file before any mutation.
    for name, digest in manifest['files'].items():
        source = PACKAGE/'pages'/name
        if Path(name).name != name or hashlib.sha256(source.read_bytes()).hexdigest() != digest:
            raise ValueError('migration package integrity mismatch: '+name)
        target = wiki_dir/name
        if name != 'Home.md' and target.exists() and (observed is None or name not in observed):
            raise ValueError('target page already exists: '+name)
    for name in manifest['preserve_existing']:
        if not (wiki_dir/name).is_file():
            raise ValueError('historical page missing: '+name)
    for name in manifest['files']:
        if apply:
            shutil.copyfile(PACKAGE/'pages'/name, wiki_dir/name)
        print(('prepared ' if apply else 'would prepare ')+name)
    return len(manifest['files'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wiki-dir', type=Path, required=True)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    prepare(args.wiki_dir, args.apply)


if __name__ == '__main__':
    main()
