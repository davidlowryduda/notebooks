# coding: utf-8
with open("bigger.dat", "r") as infile:
    xs = []
    ys = []
    for line in infile:
        x, y = line.split()
        xs.append(float(x))
        ys.append(float(y))
scatter_plot(list(zip(xs, ys)), markersize=RealNumber('0.5'), edgecolor='k', marker='.', facecolor='k')
scatter_plot(list(zip(xs, ys)), markersize=RealNumber('0.1'), edgecolor='k', marker='.', facecolor='k')
scatter_plot(list(zip(xs, ys)), markersize=RealNumber('0.1'), edgecolor='k', marker='.', facecolor='k')
scatter_plot(list(zip(xs, ys)), markersize=RealNumber('0.1'), edgecolor='k', marker='.', facecolor='k', alpha=RealNumber('0.2'))
scatter_plot(list(zip(xs, ys)), markersize=RealNumber('0.5'), edgecolor='k', marker='.', facecolor='k', alpha=RealNumber('0.2'))
scatter_plot(list(zip(xs, ys)), markersize=RealNumber('0.5'), edgecolor='k', marker='.', facecolor='k', alpha=RealNumber('0.1'))
