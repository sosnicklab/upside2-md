import numpy as np
import scipy as sp
import pandas as pd

def str_exp(x, k, b):
    """
    define stretched exponential function
    """
    return (1-np.exp(-(k*x)**b))*100 # stretched exponential

def plot_uptake(k, b, t_int):
    """
    plot D uptake curve
    """
    return (1-np.exp(-(k*t_int)**b))*100 # D uptake

def numerical_derivative(x, y):
    """
    compute the numerical derivative using finite differences
    """
    dydx = np.gradient(y, x) # calculate differences for uneven spacing
    return dydx

def normalize_se(d_norm_se):
    """
    normalize the stretched exponential data
    """
    d_norm_ses = []
    for d in d_norm_se:
        if not np.all(np.isnan(d)):
            # check whether full D Uptake was reached before final timepoint
            if np.round(np.abs(d[-1] - d[-2]), decimals=2) <= 0.01:
                # check whether full D Uptake was reached before final timepoint
                d_norms_se = np.subtract(d, d[0]) / np.subtract(np.nanmax(d), d[0])*100 #CHECKME
                d_norm_ses.append(d_norms_se)

            else:
                d_norms_se = d
                d_norm_ses.append(d_norms_se)
        else:
            d_norms_se = d
            d_norm_ses.append(d_norms_se)

    return d_norm_ses

def average_int(res_df, row_df, peps_row, peps_start_row):
    """
    calculate average integrals at residue-level
    """
    comb_res = []
    comb_res_idx = []
    comb_res_COF = []
    for index, row in res_df.iterrows():
        comb_res.extend([i for i in row['{}'.format(peps_row)]])
        comb_res_idx.extend(np.add(int(row['{}'.format(peps_start_row)]), [i for i in range(len(row['{}'.format(peps_row)]))]))
        comb_res_COF.extend([row['{}'.format(row_df)] for _ in range(len(row['{}'.format(peps_row)]))])

    return np.array(comb_res), np.array(comb_res_idx), np.array(comb_res_COF)
