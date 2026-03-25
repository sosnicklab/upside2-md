# a variety of functions and fittings
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def numerical_derivative(x, y):
    dydx = np.gradient(y, x)
    return dydx


def myexp(x, tau):
    return 1 - np.exp(-(x / tau))


def mystretchsig(x, sig, b):
    return 1 - np.exp(-np.exp(b * (x - sig) * np.log(10)))


def fitstretchsig(ys, xs, fout=False):
    if np.sum(np.isnan(ys)):
        return [-10, -10]
    xs1 = np.log10(xs)
    mydy = numerical_derivative(xs1, ys)
    sig1 = xs1[np.argmax(mydy)]
    ptab = curve_fit(
        mystretchsig,
        xs1,
        ys,
        method='trf',
        p0=[sig1, 0.5],
        bounds=([-9, 0.0], [9.0, 1.0]),
        full_output=fout,
    )
    return ptab


def pfave(s1, b1, s2, b2):
    sf = 1.7810724179874677
    return 10 ** (s2 - s1) * sf ** (1 / b1 - 1 / b2)


def relerrpfave(s1, b1, s2, b2, se1, be1, se2, be2):
    sf = 1.7810724179874677
    relb1 = be1 / (b1 ** 2)
    relb2 = be2 / (b2 ** 2)
    return np.sqrt((np.log(10)) ** 2 * (se1 ** 2 + se2 ** 2) + (np.log(sf)) ** 2 * (relb1 ** 2 + relb2 ** 2))


def dgcalc(pf, T):
    return 0.001987 * T * np.log(pf)


def dgerr(pf, pferr, T):
    return abs(0.001987 * T * pferr / pf)


def cof(beta):
    return 10 ** 4 * np.log(10) * beta / 4.0


def getkchem(sequence, pD, T):
    myT = T
    pD_corr = pD
    mol_type = "poly"

    aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    Ea_a = 14
    Ea_b = 17
    Ea_w = 19

    Ea_asp = 1
    Ea_glu = 1.083
    Ea_his = 7.5

    Kw = 15.05

    D_plus = np.power(10, -pD_corr)
    OD_minus = np.power(10, pD_corr - Kw)

    ka_poly = (np.power(10, 1.62)) / 60
    kb_poly = (np.power(10, 10.18)) / 60
    kw_poly = (np.power(10, -1.5)) / 60

    ka_oligo = ka_poly * 2.34
    kb_oligo = kb_poly * 1.35
    kw_oligo = kw_poly * 1.585
    if mol_type == "poly":
        ka = ka_poly
        kb = kb_poly
        kw = kw_poly
    elif mol_type == "oligo":
        ka = ka_oligo
        kb = kb_oligo
        kw = kw_oligo
    else:
        ka = kb = kw = -1

    dT_R = ((1 / myT) - (1 / 293)) / 0.001987
    Ft_a = np.exp(-Ea_a * dT_R)
    Ft_b = np.exp(-Ea_b * dT_R)
    Ft_w = np.exp(-Ea_w * dT_R)
    pKc_asp = -1 * np.log10(np.power(10, -4.48) * np.exp(-Ea_asp * (((1 / myT) - (1 / 278)) / 0.001987)))
    pKc_glu = -1 * np.log10(np.power(10, -4.93) * np.exp(-Ea_glu * (((1 / myT) - (1 / 278)) / 0.001987)))
    pKc_his = -1 * np.log10(np.power(10, -7.42) * np.exp(-Ea_his * (((1 / myT) - (1 / 278)) / 0.001987)))
    lambda_ma = [0, -0.54, np.log10(np.power(10, -0.90 - pD_corr) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr)) + np.power(10, 0.90 - pKc_asp) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr))), np.log10(np.power(10, -0.60 - pD_corr) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr)) + np.power(10, -0.90 - pKc_glu) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr))), -0.52, -0.22, np.log10(np.power(10, -0.80 - pD_corr) / (np.power(10, -pKc_his) + np.power(10, -pD_corr)) + np.power(10, -pKc_his) / (np.power(10, -pKc_his) + np.power(10, -pD_corr))), -0.91, -0.56, -0.57, -0.64, -0.58, 0, -0.47, -0.59, -0.437992278, -0.79, -0.739022273, -0.4, -0.41]
    rho_ma = [0, -0.46, np.log10(np.power(10, -0.12 - pD_corr) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr)) + np.power(10, 0.58 - pKc_asp) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr))), np.log10(np.power(10, -0.27 - pD_corr) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr)) + np.power(10, 0.31 - pKc_glu) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr))), -0.43, 0.218176047, np.log10(np.power(10, -0.51 - pD_corr) / (np.power(10, -pKc_his) + np.power(10, -pD_corr)) + np.power(10, -pKc_his) / (np.power(10, -pKc_his) + np.power(10, -pD_corr))), -0.59, -0.29, -0.13, -0.28, -0.13, -0.194773472, -0.27, -0.32, -0.388518935, -0.468073126, -0.3, -0.44, -0.37]
    lambda_mb = [0, 0.62, np.log10(np.power(10, 0.69 - pD_corr) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr)) + np.power(10, 0.10 - pKc_asp) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr))), np.log10(np.power(10, 0.24 - pD_corr) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr)) + np.power(10, -0.11 - pKc_glu) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr))), -0.235859464, -0.03, np.log10(np.power(10, 0.80 - pD_corr) / (np.power(10, -pKc_his) + np.power(10, -pD_corr)) + np.power(10, -0.10 - pKc_his) / (np.power(10, -pKc_his) + np.power(10, -pD_corr))), -0.73, -0.04, -0.576252728, -0.008954843, 0.49, 0, 0.06, 0.076712254, 0.37, -0.06625798, -0.701934483, -0.41, -0.27]
    rho_mb = [0, 0.55, np.log10(np.power(10, 0.60 - pD_corr) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr)) + np.power(10, -0.18 - pKc_asp) / (np.power(10, -pKc_asp) + np.power(10, -pD_corr))), np.log10(np.power(10, 0.39 - pD_corr) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr)) + np.power(10, -0.15 - pKc_glu) / (np.power(10, -pKc_glu) + np.power(10, -pD_corr))), 0.063131587, 0.17, np.log10(np.power(10, 0.83 - pD_corr) / (np.power(10, -pKc_his) + np.power(10, -pD_corr)) + np.power(10, 0.14 - pKc_his) / (np.power(10, -pKc_his) + np.power(10, -pD_corr))), -0.23, 0.12, -0.21, 0.11, 0.32, -0.24, 0.2, 0.22, 0.299550286, 0.2, -0.14, -0.11, 0.05]

    p = {"aa": aa, "lambda_ma": lambda_ma, "rho_ma": rho_ma, "lambda_mb": lambda_mb, "rho_mb": rho_mb}
    p_df = pd.DataFrame(p)

    Fa_a = []
    Fb_b = []
    Fb_w = []
    for i, s in enumerate(sequence):
        if s == 'P':
            Fa = 0
            Fb = 0
        else:
            if i == 0:
                Fa = 0
                Fb = 0
            elif i == 1:
                Fa = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_ma'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_ma'].iloc[0] - 1.32)
                Fb = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_mb'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_mb'].iloc[0] + 1.62)
            elif i == len(sequence) - 1:
                Fa = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_ma'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_ma'].iloc[0] + 0.95)
                Fb = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_mb'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_mb'].iloc[0] - 1.80)
            else:
                Fa = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_ma'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_ma'].iloc[0])
                Fb = np.power(10, p_df.loc[p_df['aa'] == sequence[i], 'lambda_mb'].iloc[0] + p_df.loc[p_df['aa'] == sequence[i - 1], 'rho_mb'].iloc[0])
        Fa_as = Fa * D_plus * ka * Ft_a
        Fa_a.append(Fa_as)

        Fb_bs = Fb * OD_minus * kb * Ft_b
        Fb_b.append(Fb_bs)

        Fb_ws = Fb * kw * Ft_w
        Fb_w.append(Fb_ws)

    Fa_a = np.array(Fa_a)
    Fb_b = np.array(Fb_b)
    Fb_w = np.array(Fb_w)

    k_chem = Fa_a + Fb_b + Fb_w
    return k_chem
