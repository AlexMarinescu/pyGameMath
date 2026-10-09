"""Build-only MkDocs integration; never imported by gem's runtime.

Virtual pages preserve repository sources. Root-relative assets are copied once;
source-code links go to GitHub. Arithmatex TeX is converted to offline MathML.
"""
import html
import posixpath
import re
from pathlib import Path
from urllib.parse import quote, urlsplit
from mkdocs.structure.files import File
from latex2mathml.converter import convert

ROOT = Path(__file__).resolve().parents[2]
REPO = 'https://github.com/AlexMarinescu/pyGameMath/blob/master/'
VIRTUAL = {
    'development/canonical-roadmap.md': 'ROADMAP.md',
    'development/benchmarks.md': 'benchmarks/README.md',
    'development/performance-policy.md': 'benchmarks/REGRESSION_POLICY.md',
    'examples/hdr-reference.md': 'examples/hdr_sh/README.md',
    'examples/regeneration.md': 'examples/showcase/README.md',
}
ASSETS = {}
# Same syntax used by the existing source documentation's link checker.
LINK = re.compile(r'(!?\[[^\]\n]*\]\()([^\s)]+)(\))')


def on_files(files, config):
    ASSETS.clear()
    for dest, source in VIRTUAL.items():
        files.append(File.generated(config, dest, abs_src_path=str(ROOT/source)))
    # Only referenced repository images are added. No entire source-tree copying.
    for file in list(files):
        if not file.is_documentation_page() or file.inclusion.is_excluded():
            continue
        origin = ROOT/VIRTUAL[file.src_uri] if file.src_uri in VIRTUAL else Path(config.docs_dir)/file.src_uri
        for match in LINK.finditer(file.content_string):
            if not match.group(1).startswith('!'):
                continue
            uri = urlsplit(match.group(2))
            if uri.scheme or uri.netloc:
                raise ValueError('remote image dependency: '+match.group(2))
            source = (origin.parent/uri.path).resolve()
            if source.is_relative_to(ROOT/'docs'):
                continue
            if not source.is_relative_to(ROOT) or not source.is_file():
                raise ValueError('invalid repository image: '+str(source))
            dest = 'assets/repository/'+source.relative_to(ROOT).as_posix()
            if source not in ASSETS:
                ASSETS[source] = dest
                files.append(File.generated(config, dest, abs_src_path=str(source)))
    return files


def on_page_markdown(markdown, page, config, files):
    origin = ROOT/VIRTUAL[page.file.src_uri] if page.file.src_uri in VIRTUAL else Path(config.docs_dir)/page.file.src_uri
    def link(match):
        target = match.group(2)
        uri = urlsplit(target)
        if uri.scheme or uri.netloc or not uri.path:
            return match.group(0)
        source = (origin.parent/uri.path).resolve()
        if not source.is_relative_to(ROOT):
            raise ValueError('repository link escape: '+target)
        suffix = ('?'+uri.query if uri.query else '')+('#'+uri.fragment if uri.fragment else '')
        if source in ASSETS:
            dest = ASSETS[source]
        elif source in {ROOT/value for value in VIRTUAL.values()}:
            dest = next(key for key, value in VIRTUAL.items() if ROOT/value == source)
        elif source == ROOT/'docs/README.md':
            dest = 'index.md'
        elif source == ROOT/'ROADMAP.md':
            dest = 'development/canonical-roadmap.md'
        elif source.is_relative_to(ROOT/'docs'):
            dest = source.relative_to(ROOT/'docs').as_posix()
            if files.get_file_from_path(dest) is None or dest.startswith('wiki-migration/pages/'):
                return match.group(1)+REPO+quote(source.relative_to(ROOT).as_posix())+suffix+match.group(3)
        else:
            if not source.exists():
                raise ValueError('missing repository source: '+target)
            return match.group(1)+REPO+quote(source.relative_to(ROOT).as_posix())+suffix+match.group(3)
        relative = posixpath.relpath(dest, posixpath.dirname(page.file.src_uri) or '.')
        return match.group(1)+relative+suffix+match.group(3)
    # Do not reinterpret links inside fenced examples or shell command strings.
    pieces = re.split(r'(```.*?```)', markdown, flags=re.DOTALL)
    return ''.join(piece if piece.startswith('```') else LINK.sub(link, piece) for piece in pieces)


def on_page_content(content, page, config, files):
    pattern = r'<(span|div) class="arithmatex">(.*?)</\1>'
    def mathml(match):
        tag, value = match.groups()
        tex = html.unescape(value)
        if tag == 'span':
            assert tex.startswith(r'\(') and tex.endswith(r'\)')
        else:
            assert tex.startswith(r'\[') and tex.endswith(r'\]')
        return '<'+tag+' class="equation">'+convert(tex[2:-2], display='block' if tag == 'div' else 'inline')+'</'+tag+'>'
    return re.sub(pattern, mathml, content, flags=re.DOTALL)
