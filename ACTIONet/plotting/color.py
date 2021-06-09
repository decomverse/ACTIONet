import colorsys
from plotly.colors import unlabel_rgb
import numpy as np

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


def hex_to_rgba(color):
    color = color.strip('#')
    # Add alpha if not included
    if len(color) < 8:
        color += 'FF'
    return tuple(int(color[i:i+2], 16)  for i in range(0, len(color), 2))


def hex_to_rgb(color):
    color = color.strip('#')
    return tuple(int(color[i:i+2], 16)  for i in range(0, len(color), 2))


def adjust_lightness(rgb, amount):
    '''
    return tuple (r,g,b) 
    '''
    if type(rgb)==str:
        rgb=unlabel_rgb(rgb)
    r, g, b = rgb
    if max(rgb)>1: 
        r=float(r)/255
        g=float(g)/255
        b=float(b)/255
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    new_r, new_g, new_b = colorsys.hls_to_rgb(h, max(0., min(1., l * amount)), s)
    r=max(0., min(1., new_r))
    g=max(0., min(1., new_g))
    b=max(0., min(1., new_b))
    return (int(r*255), int(g*255), int(b*255))


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def append_alpha_to_rgb(colors, alpha, unzip_colors=False):
    if unzip_colors:
        colors = np.array(list(zip(*colors)))
    alpha = np.atleast_2d(np.array(alpha))
    colors_rgba = np.vstack((colors, alpha)).transpose()
    colors_rgba = [tuple(x) for x in colors_rgba]
    return colors_rgba
