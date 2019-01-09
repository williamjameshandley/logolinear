import numpy
import matplotlib.pyplot as plt

def equalise(ax):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_aspect((xmax-xmin)/(ymax-ymin))

def contourplot(ax, X, Y, Z, cmap=plt.cm.cubehelix_r, *args, **kwargs):
    colorscheme = ax.contourf(X, Y, Z, cmap=cmap, *args, **kwargs) 
    for c in colorscheme.collections:
        c.set_edgecolor("face")
    cbar = ax.figure.colorbar(colorscheme, fraction=0.0455, pad=0.04)
    equalise(ax)
    return cbar

def ax_rotation(ax, p1, p2):
    tp1 = ax.transData.transform_point(p1)
    tp2 = ax.transData.transform_point(p2)
    rise = tp2[1]-tp1[1]
    run = tp2[0]-tp1[0]
    return numpy.degrees(numpy.arctan2(rise,run))


def text_along_line(ax, text, f, x, y, **kwargs):
    i = int(f * len(x))
    p1 = (x[i], y[i])
    p2 = (x[i+1], y[i+1])
    annotation = ax.annotate(text, xy=p1,verticalalignment='bottom', horizontalalignment='center', rotation_mode='anchor', **kwargs)
    annotation.set_rotation(ax_rotation(ax,p1,p2))
    return annotation

