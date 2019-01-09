import matplotlib.pyplot as plt


plt.rcParams['font.size'] = 8
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "serif"
plt.rcParams['font.serif'] = "cm"

def set_width(typ,scale=None):
    column_width = 246.0
    text_width = 510.0
    inches_per_pt = 1.0/72.27

    if typ == "text":
        fig_width = text_width * inches_per_pt 
        if scale is None:
            scale = 0.7
        fig_height = fig_width * scale

    elif typ == "column":
        fig_width = column_width * inches_per_pt 
        if scale is None:
            scale = 1.0
        fig_height = fig_width * scale

    else:
        raise ValueError("Unknown width, should be \"text\" or \"column\"")

    plt.rcParams['figure.figsize'] = [fig_width, fig_height]

set_width("column")
