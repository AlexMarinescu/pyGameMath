"""Validate emitted documentation, not only source Markdown (standard library)."""
import argparse
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.ids, self.links, self.images, self.resources = set(), [], [], []
        self.math, self.highlighted, self.nav = 0, False, []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if 'id' in a:
            self.ids.add(a['id'])
        if 'href' in a:
            self.links.append(a['href'])
            if 'md-nav__link' in a.get('class', ''):
                self.nav.append(a['href'])
            if tag == 'link' and a.get('rel') in ('stylesheet', 'icon'):
                self.resources.append(a['href'])
        if 'src' in a:
            self.links.append(a['src']); self.resources.append(a['src'])
            if tag == 'img':
                self.images.append(a)
        if tag == 'math':
            self.math += 1
        if tag == 'span' and a.get('class') in ('kn', 'k', 'nf'):
            self.highlighted = True


def check(site_dir):
    root = Path(site_dir).resolve()
    pages = {p.resolve(): Page(p.read_text()) for p in root.rglob('*.html')}
    if not pages:
        raise ValueError('no built HTML pages')
    count = images = 0
    for path, page in pages.items():
        for target in page.links:
            uri = urlsplit(target)
            if uri.scheme or uri.netloc:
                continue
            if uri.path.startswith('/'):
                raise ValueError('hosting-prefix-unsafe absolute URL: '+target)
            dest = (path.parent/unquote(uri.path)).resolve() if uri.path else path
            if dest.is_dir():
                dest /= 'index.html'
            if not dest.is_relative_to(root) or not dest.is_file():
                raise ValueError(str(path.relative_to(root))+': missing link '+target)
            if uri.fragment and dest in pages and unquote(uri.fragment) not in pages[dest].ids:
                raise ValueError(str(path.relative_to(root))+': missing fragment '+target)
            count += 1
        for asset in page.resources:
            uri = urlsplit(asset)
            if uri.scheme not in ('', 'data') or uri.netloc:
                raise ValueError('external asset dependency: '+asset)
        for image in page.images:
            if not image.get('alt'):
                raise ValueError('missing image alternative text: '+str(path))
            images += 1
    for css in root.rglob('*.css'):
        for target in re.findall(r'url\([\"\']?([^\)\"\']+)', css.read_text()):
            uri = urlsplit(target)
            if uri.scheme == 'data':
                continue
            if uri.scheme or uri.netloc or uri.path.startswith('/'):
                raise ValueError('non-local CSS asset: '+target)
            if not (css.parent/unquote(uri.path)).is_file():
                raise ValueError('missing CSS asset: '+target)
    required = ['api/vector/index.html', 'api/matrix/index.html', 'api/quaternion/index.html',
                'api/spherical-harmonics/index.html', 'examples/index.html',
                'development/canonical-roadmap/index.html', 'architecture/notation/index.html']
    for name in required:
        if root/name not in pages:
            raise ValueError('missing required page: '+name)
    nav = pages[root/'index.html'].nav
    for name in required:
        if name[:-10] not in nav and name.replace('index.html', '') not in nav:
            raise ValueError('required page absent from homepage navigation: '+name)
    search = json.loads((root/'search/search_index.json').read_text())
    api = [d for d in search['docs'] if d['location'].startswith('api/')]
    corpus = '\n'.join(d['text'] for d in api)
    for symbol in ['rotate_coefficients', 'quat_slerp', 'findDrawingPoints', 'viewport']:
        if symbol not in corpus:
            raise ValueError('API missing from search index: '+symbol)
    notation = pages[root/'architecture/notation/index.html']
    if notation.math < 10 or not notation.highlighted:
        raise ValueError('equations or Python syntax highlighting absent')
    return dict(html_pages=len(pages), local_links_checked=count, image_references=images,
                search_api_entries=len(api), mathml_equations=notation.math,
                syntax_highlighting=True, external_asset_dependencies=0,
                navigation_checked=len(required), hosting='relative project-prefix-safe paths')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--site-dir', type=Path, default=Path('site'))
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    result = check(a.site_dir)
    if a.output:
        a.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
