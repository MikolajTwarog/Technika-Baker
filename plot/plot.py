import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import numpy
import sys

ARGS = len(sys.argv) - 2
IN_FILE1 = 'results/' + sys.argv[1]
IN_FILE2 = 'results/' + sys.argv[2]
OUT_FILE = 'results/' + sys.argv[-2]
X_AX = sys.argv[-1]
Y_AX = 'czas w ms'

colors = ['#329FBB', '#D1431F', '#8AC272']
linestyles = ['-', '--', '-.', ':']

def figure_size(figure_size_scale):
    inches_per_pt = 1.0 / 72.27 # Convert pt to inch
    golden_mean = (numpy.sqrt(5.0) - 1.0) / 2.0 # Aesthetic ratio (you could change this)
    figure_width = 455.24411 * inches_per_pt * figure_size_scale
    figure_height = figure_width * golden_mean
    return [figure_width, figure_height]

publication_with_latex = {
        "pgf.texsystem": "pdflatex", # change this if using xetex or lautex
        # "text.usetex": True, # use LaTeX to write all text
        "font.family": "serif",
        "axes.labelsize": 8, # LaTeX default is 10pt font.
        "font.size": 8,
        "legend.fontsize": 8, # Make the legend/label fonts a little smaller
        "savefig.dpi": 125,
        "text.latex.preamble": r"\usepackage{amsmath,amssymb,amsfonts}",
        "figure.figsize": figure_size(1)
    }
plt.style.use('ggplot')
matplotlib.rcParams.update(publication_with_latex)

fig, ax = plt.subplots()
big = False

for i in range(1, ARGS):
    T1 = []
    T2 = []
    with open('results/' + sys.argv[i], "r") as f:
        name = f.readline()
        for line in f.readlines():
            T1.append(int(line.split(" ")[0]))
            T2.append(float(line.split(" ")[1]))
    if (len(T1) > 15):
        T3 = numpy.polyfit(T1, T2, 1)
        p = numpy.poly1d(T3)
        xp = numpy.linspace(0, T1[-1], 100)
        ax.plot(xp, p(xp), label=name, color=colors[i-1], linestyle=linestyles[i-1])
        ax.scatter(T1, T2, color=colors[i-1], alpha=0.3, s=5)
        big = True
    else:
        ax.plot(T1, T2, label=name, color=colors[i-1], linestyle=linestyles[i-1])

if big:
    ax.set_ylim(bottom=1)
ax.legend(loc='lower right', frameon=True)
ax.set_yscale('log')

plt.ylabel(Y_AX)
plt.xlabel(X_AX)
plt.savefig(OUT_FILE)