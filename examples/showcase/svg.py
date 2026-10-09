"""Deterministic SVG drawing primitives (screen coordinates, CSS sRGB colors).

Only geometry is generated: generic sans-serif text uses the viewer's fonts.
All text and attributes are XML escaped; no scripts or external resources.
"""
from html import escape

INK = '#18344a'
BLUE = '#176cba'
ORANGE = '#bf4d16'
GREEN = '#137b62'
MUTED = '#667987'


def number(value):
    return format(0.0 if abs(value) < 5e-10 else value, '.6f').rstrip('0').rstrip('.')


class Canvas:
    def __init__(self, title, subtitle):
        self.parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="960" height="560" viewBox="0 0 960 560" role="img">',
                      '<title>'+escape(title)+'</title>', '<desc>'+escape(subtitle)+'</desc>',
                      '<rect width="960" height="560" fill="#f5f8fb"/>']
        self.text(32, 40, title, 25)
        self.text(32, 69, subtitle, 15, MUTED)

    def text(self, x, y, text, size=15, color=INK):
        self.parts.append('<text x="{}" y="{}" fill="{}" font-family="sans-serif" font-size="{}">{}</text>'.format(number(x), number(y), escape(color, quote=True), size, escape(str(text))))

    def line(self, a, b, color=MUTED, width=2, dashed=False):
        self.parts.append('<path d="M {} {} L {} {}" fill="none" stroke="{}" stroke-width="{}"{}/>'.format(number(a[0]), number(a[1]), number(b[0]), number(b[1]), escape(color, quote=True), width, ' stroke-dasharray="6 5"' if dashed else ''))

    def arrow(self, a, b, color=BLUE):
        import math
        self.line(a, b, color, 3)
        angle = math.atan2(b[1]-a[1], b[0]-a[0])
        for offset in (-.5, .5):
            self.line(b, (b[0]-11*math.cos(angle+offset), b[1]-11*math.sin(angle+offset)), color, 3)

    def circle(self, point, radius=4, color=BLUE):
        self.parts.append('<circle cx="{}" cy="{}" r="{}" fill="{}"/>'.format(number(point[0]), number(point[1]), radius, escape(color, quote=True)))

    def path(self, points, color=BLUE, width=3, close=False, dashed=False):
        if not points:
            return
        data = 'M '+' L '.join(number(x)+' '+number(y) for x, y in points)
        if close:
            data += ' Z'
        self.parts.append('<path d="{}" fill="none" stroke="{}" stroke-width="{}"{}/>'.format(data, escape(color, quote=True), width, ' stroke-dasharray="6 5"' if dashed else ''))

    def finish(self):
        return '\n'.join(self.parts+['</svg>', ''])


def chart(canvas, origin, scale, extent=(-2, 5, -2, 5)):
    """Return a drawing-only XY-to-screen mapping (screen Y points down)."""
    def screen(point):
        return origin[0]+scale*point[0], origin[1]-scale*point[1]
    xmin, xmax, ymin, ymax = extent
    canvas.line(screen((xmin, 0)), screen((xmax, 0)), '#c0ccd5', 1)
    canvas.line(screen((0, ymin)), screen((0, ymax)), '#c0ccd5', 1)
    canvas.text(*screen((xmax, -.35)), 'X', 13)
    canvas.text(*screen((.15, ymax)), 'Y', 13)
    return screen
