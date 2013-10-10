import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt


mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Helvetica'
mpl.rcParams['figure.figsize'] = [8.0, 5.0]
tick_size = 15
fontlabel_size = 16
params = {'axes.labelsize': fontlabel_size, 'text.fontsize': fontlabel_size, 'legend.fontsize': fontlabel_size, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'text.usetex': True}
mpl.rcParams.update(params)

if __name__ == "__main__":
    x = np.array([1,2,3])
    y = np.array([4,5,6])
    plt.plot(x,y,'go')
    plt.xlim(0,4)
    plt.ylim(3,7)
    plt.xlabel("X Axis")
    plt.ylabel("Y Axis")
    plt.show()
