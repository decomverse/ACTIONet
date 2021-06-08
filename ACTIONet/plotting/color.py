import colorsys
from plotly.colors import unlabel_rgb 
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
