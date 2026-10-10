"""Optional local-browser QA; requires Playwright and a Chromium executable.

Not required for builds or gem imports. Starts a loopback-only static server and
also exercises a project-prefix mount. Checks that no external requests occur.
"""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
import re
from pathlib import Path
import threading
from playwright.sync_api import sync_playwright


class Handler(SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path.startswith('/pyGameMath/'):
            self.path = self.path[len('/pyGameMath'):]
        super().do_GET()

    def log_message(self, *args):
        pass


def link_contrast(page, selector='main .md-typeset p a:not(.md-button)', background='body'):
    colors = page.evaluate("([selector,background]) => [getComputedStyle(document.querySelector(selector)).color, getComputedStyle(document.querySelector(background)).backgroundColor]", [selector,background])
    def luminance(color):
        rgb = [int(v)/255 for v in re.findall(r'\d+', color)[:3]]
        linear = [v/12.92 if v<=.04045 else ((v+.055)/1.055)**2.4 for v in rgb]
        return sum(a*b for a,b in zip(linear,[.2126,.7152,.0722]))
    a,b = sorted(map(luminance,colors))
    ratio = (b+.05)/(a+.05)
    assert ratio >= 4.5, (colors,ratio)
    return ratio


def set_scheme(page, scheme):
    page.evaluate("scheme => document.body.setAttribute('data-md-color-scheme',scheme)", scheme)
    # Measure the settled palette, not a partially animated link color.
    page.wait_for_function("""() => {
      const link = document.querySelector('main .md-typeset p a:not(.md-button)');
      const probe = document.createElement('span');
      probe.style.color = getComputedStyle(document.body).getPropertyValue('--md-typeset-a-color');
      document.body.appendChild(probe);
      const expected = getComputedStyle(probe).color;
      probe.remove();
      return link && getComputedStyle(link).color === expected;
    }""")


def check(site_dir, executable, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    server = ThreadingHTTPServer(('127.0.0.1', 0), partial(Handler, directory=str(site_dir.resolve())))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    origin = 'http://127.0.0.1:'+str(server.server_port)
    external, errors, missing = [], [], []
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(executable_path=executable, headless=True)
            context = browser.new_context(viewport={'width': 1440, 'height': 1000})
            context.on('request', lambda request: external.append(request.url) if not request.url.startswith((origin+'/', 'data:', 'blob:')) else None)
            page = context.new_page()
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.on('response', lambda response: missing.append(response.url) if response.status >= 400 else None)
            page.goto(origin+'/', wait_until='networkidle')
            assert page.title().startswith('pyGameMath')
            assert page.locator('main h1').inner_text().startswith('Graphics mathematics')
            assert page.locator('main img').count() == 3
            assert page.locator('main img').evaluate_all('(images) => images.every(i => i.complete && i.naturalWidth > 0)')
            contrasts = {'light': link_contrast(page)}
            brand_contrast = {
                'py': link_contrast(page, '.md-header .gem-wordmark__py', '.md-header'),
                'GameMath': link_contrast(page, '.md-header .gem-wordmark__name', '.md-header')}
            assert page.locator('.md-header .gem-wordmark').text_content() == 'pyGameMath'
            assert page.locator('.gem-card').count() == 4
            assert page.locator('.md-tabs a').filter(has_text='Tutorials').count() == 1
            version = browser.version
            page.screenshot(path=str(output_dir/'home-light.png'), full_page=True)
            page.screenshot(path=str(output_dir/'home-light-viewport.png'))
            page.locator('label[title="Switch to dark mode"]').click()
            assert page.locator('body').get_attribute('data-md-color-scheme') == 'slate'
            set_scheme(page, 'slate')
            contrasts['dark'] = link_contrast(page)
            page.screenshot(path=str(output_dir/'home-dark.png'), full_page=True)
            page.screenshot(path=str(output_dir/'home-dark-viewport.png'))
            page.goto(origin+'/architecture/notation/', wait_until='networkidle')
            assert page.locator('math').count() == 24
            assert page.locator('math mtable').count() >= 1
            assert page.locator('math').evaluate_all('(items) => items.every(i => i.getBoundingClientRect().height > 0)')
            page.screenshot(path=str(output_dir/'notation-dark.png'), full_page=True)
            # Search through the actual theme UI, not only the index JSON.
            search = page.get_by_role('textbox', name='Search')
            search.click()
            search.press_sequentially('rotate_coefficients', delay=40)
            try:
                page.wait_for_function("document.querySelector('.md-search-result__list')?.textContent.includes('rotate_coefficients')", timeout=10000)
            except Exception:
                print('Search diagnostic:', page.locator('.md-search-result__meta').inner_text(), errors, missing, external)
                page.screenshot(path=str(output_dir/'search-diagnostic.png'))
                raise
            assert page.locator('.md-search-result__list a[href*="spherical-harmonics"]').count() > 0
            page.keyboard.press('Escape')
            page.goto(origin+'/examples/', wait_until='networkidle')
            assert page.locator('main img').evaluate_all('(images) => images.every(i => i.naturalWidth > 0 && Math.abs(i.width/i.height-i.naturalWidth/i.naturalHeight)<.01)')
            page.screenshot(path=str(output_dir/'gallery-dark.png'), full_page=True)
            page.set_viewport_size({'width': 390, 'height': 844})
            page.goto(origin+'/pyGameMath/', wait_until='networkidle')
            assert page.evaluate('document.documentElement.scrollWidth <= innerWidth+1')
            page.screenshot(path=str(output_dir/'home-mobile.png'), full_page=True)
            menu = page.locator('.md-header__button[for="__drawer"]')
            menu.click()
            assert page.locator('#__drawer').is_checked()
            page.locator('.md-nav--primary').hover()
            page.mouse.wheel(0, -10000)
            page.screenshot(path=str(output_dir/'mobile-drawer.png'))
            page.locator('.md-nav--primary .md-nav__container').filter(
                has=page.get_by_role('link',name='Getting started',exact=True)).locator('label').click()
            assert page.locator('.md-nav--primary a[href*="getting-started/installation"]').count() == 1
            page.locator('.md-nav--primary a[href*="getting-started/installation"]').click()
            page.wait_for_url('**/pyGameMath/getting-started/installation/')
            page.wait_for_load_state('networkidle')
            assert '/pyGameMath/getting-started/installation/' in page.url
            assert page.locator('main h1').inner_text().startswith('Install current development code')
            assert not errors and not missing and not external, (errors, missing, external)
            # Basic semantic accessibility checks; not a full WCAG audit.
            assert page.locator('html').get_attribute('lang') == 'en'
            assert page.locator('main').count() == 1
            assert page.locator('a.md-skip').count() == 1
            # Tables and equations may scroll locally, never the whole viewport.
            layouts = []
            for scheme in ('default', 'slate'):
                set_scheme(page, scheme)
                for width in (320, 390, 768, 1440):
                    page.set_viewport_size({'width': width, 'height': 900})
                    for route in ('', 'api/vector/', 'api/quaternion/', 'api/legendre/',
                                  'architecture/notation/', 'tutorials/bezier/', 'examples/'):
                        page.goto(origin+'/pyGameMath/'+route, wait_until='networkidle')
                        # Keep the selected scheme across each full-page visit.
                        set_scheme(page, scheme)
                        assert page.evaluate('document.documentElement.scrollWidth <= innerWidth+1'), (scheme,width,route)
                        assert page.locator('math').evaluate_all('(items) => items.every(i => i.getBoundingClientRect().height > 0)')
                        assert page.locator('main img').evaluate_all('(items) => items.every(i => i.complete && i.naturalWidth > 0)')
                        layouts.append({'scheme':scheme,'width':width,'route':route,'overflow':False})
                    page.goto(origin+'/pyGameMath/', wait_until='networkidle')
                    set_scheme(page, scheme)
                    assert link_contrast(page) >= 4.5
                    assert link_contrast(page, 'main .md-typeset p') >= 4.5
                    assert link_contrast(page, 'main figcaption') >= 4.5
                    assert link_contrast(page, 'main .md-button--primary', 'main .md-button--primary') >= 4.5
            page.set_viewport_size({'width':1440,'height':1000})
            page.goto(origin+'/',wait_until='networkidle')
            page.locator('a.md-skip').focus()
            assert page.locator('a.md-skip').evaluate('(node) => node === document.activeElement')
            page.keyboard.press('Enter')
            assert page.url.endswith('#graphics-mathematics-in-python')
            assert not errors and not missing and not external, (errors, missing, external)
            browser.close()
    finally:
        server.shutdown(); server.server_close()
    return dict(browser='Chromium '+version, link_contrast_ratios=contrasts, desktop=[1440,1000], mobile=[390,844],
                light_dark=True, search='rotate_coefficients result verified',
                mathml=True, mathml_notation_equations=24, gallery_aspect_ratios=True, mobile_drawer=True,
                brand_contrast_ratios=brand_contrast, layouts_checked=layouts, keyboard_skip_link=True,
                project_prefix='/pyGameMath/', external_requests=external,
                javascript_errors=errors, missing_resources=missing,
                accessibility='language, main landmark, skip link, controls and alt text; not a full WCAG audit',
                screenshots=[p.name for p in sorted(output_dir.glob('*.png'))])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site-dir', type=Path, default=Path('site'))
    parser.add_argument('--browser-executable', default='/usr/bin/chromium')
    parser.add_argument('--output-dir', type=Path, default=Path('/tmp/gem-site-browser'))
    args = parser.parse_args()
    result = check(args.site_dir, args.browser_executable, args.output_dir)
    (args.output_dir/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
