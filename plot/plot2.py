import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy
import sys


OUT_FILE = 'results/is_ptas_res.pdf'
X_AX = 'k'
Y_AX = 'rozmiar zbioru niezależnego'

colors = ['#329FBB', '#D1431F', '#8AC272']

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

IN_FILE1 = sys.argv[1]
T1 = []
T2 = []
T3 = []
T4 = []
with open(IN_FILE1, "r") as f:
    for line in f.readlines():
        T1.append(int(line.split(" ")[0]))
        T2.append(float(line.split(" ")[1]))
        T3.append(int(line.split(" ")[0]))
        T4.append((float(line.split(" ")[1]) * (int(line.split(" ")[0]) + 1)) / int(line.split(" ")[0]))

fig, ax = plt.subplots()
l1 = ax.bar(T3, T4, label="maksymalne rozwiązanie optymalne", hatch="/////", edgecolor='lightblue', color=colors[0])
l2 = ax.bar(T1, T2, label="rozwiązanie techniki Baker", color=colors[1])
# ax.fill_between(T1, T2, T4, color='lightblue')
legend = [Patch(edgecolor='lightblue', label='Color Patch'), Patch()]
ax.legend((l1, l2), ("maksymalne rozwiązanie optymalne", "rozwiązanie techniki Baker"), loc='lower right', frameon=True)

plt.ylabel(Y_AX)
plt.xlabel(X_AX)
plt.savefig(OUT_FILE)