import colorsys

def hex_to_rgba(color):
    color = color.strip('#')
    # Add alpha if not included
    if len(color) < 8:
        color += 'FF'
    return tuple(int(color[i:i+2], 16) / 255 for i in range(0, len(color), 2))

def adjust_lightness(rgba, amount):
    r, g, b, a = rgba
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    new_r, new_g, new_b = colorsys.hls_to_rgb(h, max(0., min(1., l * amount)), s)
    return (max(0., min(1., new_r)), max(0., min(1., new_g)), max(0., min(1., new_b)), a)
