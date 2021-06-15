import matplotlib
import matplotlib.pyplot as plt
import numpy
import sys


IN_FILE1 = 'results/' + sys.argv[1]
IN_FILE2 = 'results/' + sys.argv[2]
OUT_FILE = 'results/' + sys.argv[3]
X_AX = sys.argv[4]
Y_AX = 'czas w ms'

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

T1 = []
T2 = []
T3 = []
T4 = []
with open(IN_FILE1, "r") as f:
    for line in f.readlines():
        T1.append(int(line.split(" ")[0]))
        T2.append(float(line.split(" ")[1]))
        
with open(IN_FILE2, "r") as f:
    for line in f.readlines():
        T3.append(int(line.split(" ")[0]))
        T4.append(float(line.split(" ")[1]))

fig, ax = plt.subplots()
ax.plot(T3, T4, label="Bodlaender")
ax.plot(T1, T2, label="Baker")
ax.legend(loc='upper left', frameon=False)
ax.set_yscale('log')

plt.ylabel(Y_AX)
plt.xlabel(X_AX)
plt.savefig(OUT_FILE)