import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.colors as colors
from matplotlib.patches import Circle
import seaborn as sns
from scipy.stats import pearsonr
import os

plt.rcParams['mathtext.fontset'] = "stix"

# MODELS
def linear(x, params, pivots):
    b = params[0]
    m = params[1]
    xp = pivots[0]
    return b + m*(x - xp)

def dLinear(x, params, pivots):
    b = params[0]
    m = params[1]

    return m

def bhc2(x, params, pivots):
    b1BH = params[0]
    theta1BH = params[1]

    b2BH = params[2]
    theta2BH = params[3]

    pivot1 = pivots[0]
    pivot2 = pivots[1]

    c2 = x

    return -np.log(np.exp(-b1BH - np.tan(theta1BH)*(c2 - pivot1)) + np.exp(-b2BH - np.tan(theta2BH)*(c2 - pivot2)))

def dbhc2(x, params, pivots):
    b1BH = params[0]
    theta1BH = params[1]

    b2BH = params[2]
    theta2BH = params[3]

    pivot1 = pivots[0]
    pivot2 = pivots[1]

    c2 = x

    top = -np.exp(-b1BH - np.tan(theta1BH)*(c2 - pivot1))*np.tan(theta1BH) - np.exp(-b2BH - np.tan(theta2BH)*(c2 - pivot2))*np.tan(theta2BH)
    bottom = np.exp(-b1BH - np.tan(theta1BH)*(c2 - pivot1)) + np.exp(-b2BH - np.tan(theta2BH)*(c2 - pivot2))

    return -top/bottom

def rvc2(x, params, pivots):
    b1RV = params[0]
    theta1RV = params[1]

    b2RV = params[2]
    theta2RV = params[3]

    pivot1 = pivots[0]
    pivot2 = pivots[1]

    c2 = x

    return np.log(np.exp(b1RV + np.tan(theta1RV)*(c2 - pivot1)) + np.exp(b2RV + np.tan(theta2RV)*(c2 - pivot2)))

def drvc2(x, params, pivots):
    b1RV = params[0]
    theta1RV = params[1]

    b2RV = params[2]
    theta2RV = params[3]

    pivot1 = pivots[0]
    pivot2 = pivots[1]

    c2 = x

    top = np.exp(b1RV + np.tan(theta1RV)*(c2 - pivot1))*np.tan(theta1RV) + np.exp(b2RV + np.tan(theta2RV)*(c2 - pivot2))*np.tan(theta2RV)
    bottom = np.exp(b1RV + np.tan(theta1RV)*(c2 - pivot1)) + np.exp(b2RV + np.tan(theta2RV)*(c2 - pivot2))

    return top / bottom

# PARAMETERS

hist_filenames = {
    "testlin":"testlin_sloppy_pivot_1000000.txt",
    "c1c2":"c1c2_100000.txt",
    "bhc2":"bhc2_100000.txt",
    "rvc2":""
}

data_filenames = {
    "testlin":"testlin_sloppy_pivot_1000000.txt",
    "c1c2":"c1c2_data.csv",
    "bhc2":"bhc2_data.csv",
    "rvc2":"rvc2_data.csv"
}

xlabels = {
    "testlin":"$x$",
    "c1c2":"$c_2$",
    "bhc2":"$c_2$",
    "rvc2":"$c_2$"
}

ylabels = {
    "testlin":"$y$",
    "c1c2":"$c_1$",
    "bhc2":"BH",
    "rvc2":"$R_V$"
}

best_fit_params = { # including slop
    "testlin":[0.870858, 0.978703, 0.602990, 0.585983],
    "c1c2":[-0.905186, -3.151891, 0.059483, 0.171067],
    "bhc2":[1.983692, 4.567873, 2.403995, -1.134415, 0.129025, 0.495434],
    "rvc2":"$R_V$"
}

all_pivots = {
    # "testlin":
    "c1c2":[1.024747],
    "bhc2":[-0.086661, 1.334268],
    # "rvc2":"$R_V$"
}

all_param_names = {
    "testlin":["$b$", "$m$", "$\sigma_x$", "$\sigma_y$"],
    "c1c2":["$b^{c_1}$", "$m^{c_1}$", "$\sigma_{c_2}^{c_1}$", "$\sigma_{c_1}$"],
    "bhc2":["$b_1^{BH}$", "$\\theta_1^{BH}$", "$b_2^{BH}$", "$\\theta_2^{BH}$", "$\sigma_{c_2}^{BH}$", "$\sigma_{BH}$"],
    "rvc2":["$b_1^{R_V}$", "$\\theta_1^{R_V}$", "$b_2^{R_V}$", "$\\theta_2^{R_V}$", "$\sigma_{c_2}^{R_V}$", "$\sigma_{R_V}$"]
}

models = {
    "testlin":linear,
    "c1c2":linear,
    "bhc2":bhc2,
    "rvc2":rvc2
}

dModels = {
    "testlin":dLinear,
    "c1c2":dLinear,
    "bhc2":dbhc2,
    "rvc2":drvc2
}

all_param_plot_ranges_1D = {
    "testlin":((0,2), (4,8)),
    "c1c2":((-1.1, -0.8), (-3.4, -2.8), (0.0525, 0.065), (0.065, 0.265)),
    "bhc2":((1.96, 2.0), (4.565, 4.58), (2.404, 2.410), (-1.25, -1.1), (0.115, 0.14), (0.48, 0.53)),
    # "rvc2":"$R_V$"
}

all_param_plot_ranges_2D = all_param_plot_ranges_1D
# all_param_plot_ranges_2D = {
#     "testlin":((0,2), (4,8)),
#     "c1c2":((-1.1, -0.7), (-3.3, -3),  (0.05, 0.07), (0.16, 0.18))
#     # ,"bhc2":"BH",
#     # "rvc2":"$R_V$"
# }

data_plot_ranges = {
    #"testlin":((0,2), (4,8)),
    "c1c2":((-1, 3), (-8, 7)),
    "bhc2":((-1, 3), (-1, 8))
    # "rvc2":"$R_V$"
}



if __name__ == '__main__':
    # print("Make sure to translate slope to theta")

    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("modelname", type=str, help="model name")
    args = parser.parse_args()
    modelname = args.modelname

    filename = os.path.join("./histdata/",hist_filenames[modelname])
    datafilename = os.path.join("./rawdata/",data_filenames[modelname])

    xlabel = xlabels[modelname]
    ylabel = ylabels[modelname]

    fitted_params = best_fit_params[modelname]
    model = models[modelname]
    dModel = dModels[modelname]
    param_names = all_param_names[modelname]
    pivots = all_pivots[modelname]

    param_plot_ranges_1D = all_param_plot_ranges_1D[modelname]
    param_plot_ranges_2D = all_param_plot_ranges_2D[modelname]
    data_plot_range = data_plot_ranges[modelname]

    M = len(fitted_params)

    datadf = pd.read_csv(datafilename, sep=",", names=['x', 'sx', 'y', 'sy', 'w'])
    df = pd.read_csv(filename, sep=" ", header=None)

    # print(datadf)

    x = datadf['x']
    sx = datadf['sx']
    y = datadf['y']
    sy = datadf['sy']

    ar = df.values
    ar = np.transpose(ar)

    param_samples = []
    for i in range(M):
        param_samples.append(ar[i])

    # a0new = []
    # a1new = []
    # for i, _ in enumerate(a0):
    #     if a0[i] >= a0rng[0] and a0[i] <= a0rng[1] and a1[i] >= a1rng[0] and a1[i] <= a1rng[1]:
    #         a0new.append(a0[i])
    #         a1new.append(a1[i])
    
    # a0 = np.array(a0new)
    # a1 = np.array(a1new)

    # print(pearsonr(a0, a1))
    

    R = param_samples[0].size
    bincount = np.int(np.sqrt(float(R)))

    datawidth = x.ptp()
    x_m = np.linspace(np.min(x) - datawidth/3, np.max(x) + datawidth/3, 1000)
    y_m = model(x_m, fitted_params, pivots)

    shiftY = np.sqrt(np.power(fitted_params[-1], 2.0) + np.power(dModel(x_m, fitted_params, pivots) * fitted_params[-2], 2.0))

    plt.rc('axes', labelsize=12)     # fontsize of the axes label

    # plt.figure(num=1,figsize=(8,6),dpi=100,facecolor='white')
    plt.figure(num=1,figsize=(8,8),dpi=100,facecolor='white')
    # plot single-parameter histograms
    inds = [1,3,4,6,7,9]
    for j in range(M):
        plt.subplot(M/2, 3, inds[j])
        n, bins, patches = plt.hist(param_samples[j], bincount, density=True, facecolor='b', alpha=0.5)
        plt.xlabel(param_names[j])
        plt.ylabel("PDF")
        plt.xlim(param_plot_ranges_1D[j])
        plt.xticks(rotation=45)
        plt.yticks(rotation=45)

    # 2-param histograms
    for j in range(M//2):
        plt.subplot(M/2, 3, 2 + 3*j)
        # Kernel density/contour plot
        sns.kdeplot(param_samples[2*j], param_samples[2*j + 1], n_levels=3, shade=True, shade_lowest=False)
        # circle = plt.Circle((bfit, mfit), 0.05, color='r')
        # ax.add_artist(circle)
        plt.xlabel(param_names[2*j])
        plt.ylabel(param_names[2*j + 1])
        plt.xlim(param_plot_ranges_2D[2*j])
        plt.ylim(param_plot_ranges_2D[2*j + 1])
        plt.xticks(rotation=45)
        plt.yticks(rotation=45)

    plt.tight_layout()
    plt.show()

    # FONT SIZES
    SMALL_SIZE = 12
    MEDIUM_SIZE = 16
    BIGGER_SIZE = 18

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    # plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    # plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    #

    # plot in data space
    plt.figure(num=1,figsize=(8,8),dpi=100,facecolor='white')

    plt.errorbar(x, y, yerr=sy, xerr=sx,  fmt="o", color="k", ecolor='k', label="Data", markersize=1, elinewidth=1, alpha=0.5)
    plt.plot(x_m, y_m, label="Model", color="k", linewidth=2, alpha=0.75)
    plt.fill_between(x_m, y_m + shiftY, y_m - shiftY, color="b", alpha=0.4)
    plt.fill_between(x_m, y_m + 2*shiftY, y_m - 2*shiftY, color="b", alpha=0.3)
    plt.fill_between(x_m, y_m + 3*shiftY, y_m - 3*shiftY, color="b", alpha=0.2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(data_plot_range[0])
    plt.ylim(data_plot_range[1])
    plt.show()