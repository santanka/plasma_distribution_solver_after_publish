import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
import matplotlib as mpl
from joblib import Parallel, delayed
import multiprocessing

N_CPU = multiprocessing.cpu_count()

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern Roman']
mpl.rcParams['mathtext.fontset'] = 'cm'

plt.rcParams["font.size"] = 20

mp.mp.dps = 50

def hyp2f1(a, b, c, z):
    return mp.hyper([a, b], [c], z)

def incomplete_beta_Bz(a, b, z):
    return mp.betainc(a, b, 0, z, regularized=False)

def nV_over_Zvth_kappa(Bi, Bb, Be, kappa, T_anisotropy):
    k = mp.mpf(kappa)
    Bb = mp.mpf(Bb); Be = mp.mpf(Be); Bi = mp.mpf(Bi)
    T_anisotropy = mp.mpf(T_anisotropy)

    if Bi > Be:
        Bi = Be

    # --- 第1項（最新式：2k/(4k+1)）---
    z1 = 1. - 2. * (k - 1) / (2 * k - 1) * T_anisotropy * (Be / Bb - 1)
    term1 = mp.sqrt(k) / 2 / k * mp.gamma(k+1) / mp.gamma(k+mp.mpf('0.5')) * (Bi/Bb) * hyp2f1(k, 1, 2*k+1, z1)

    val = term1 / 2. / mp.sqrt(mp.pi)

    # 微小虚部対策
    if isinstance(val, mp.mpc):
        val = mp.re(val)

    return val

def n_over_Z_kappa_array_parallel(Bi_array, Bb, Be, kappa, T_anisotropy, n_jobs=-1):
    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(nV_over_Zvth_kappa)(Bi, Bb, Be, kappa, T_anisotropy)
        for Bi in Bi_array
    )
    return np.array([float(val) for val in results])

def nV_over_Zvth_biMaxwell(Bi, Bb, Be, T_anisotropy):
    return Bi / ((1. - T_anisotropy) * Bb + T_anisotropy * Be) / 2. / mp.sqrt(mp.pi)



# ---- parameters ----
Bb = 1.0
Be = 400.0
kappa_list = [2., 3., 5., 10., 100.]
T_anisotropy_list = [1/5, 1/2, 1., 2., 5.]  # = Tperp / Tpara

Bi_grid = np.logspace(np.log10(Bb), np.log10(Be), 100)
Bi_grid[0]  = Bb
Bi_grid[-1] = Be

# ---- compute ----
# shape: [iT][ik][iBi]
nV_kappa = [[None for _ in kappa_list] for _ in T_anisotropy_list]
nV_bimax = [None for _ in T_anisotropy_list]

for iT, T_anisotropy in enumerate(T_anisotropy_list):
    nV_bimax[iT] = nV_over_Zvth_biMaxwell(Bi_grid, Bb, Be, T_anisotropy)
    for ik, kappa in enumerate(kappa_list):
        nV_kappa[iT][ik] = n_over_Z_kappa_array_parallel(
            Bi_grid, Bb, Be, kappa, T_anisotropy, n_jobs=N_CPU
        )

# ---- plot: 1 panel per T_anisotropy, curves are kappa ----
fig = plt.figure(figsize=(10, 20), dpi=100)
color_list = ['b', 'r', 'g', 'm', 'k']

for iT, T_anisotropy in enumerate(T_anisotropy_list):
    ax = fig.add_subplot(len(T_anisotropy_list), 1, iT+1)

    for ik, kappa in enumerate(kappa_list):
        ax.plot(Bi_grid, nV_kappa[iT][ik], c=color_list[ik], lw=2, alpha=0.8,
                label=rf"$\kappa={kappa}$")

    ax.plot(Bi_grid, nV_bimax[iT], c='0.3', lw=2, ls='--', label="bi-Maxwellian")

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(left=1)
    ax.set_xlabel(r"$B ( r_{\parallel i} ) / B ( r_{\parallel b} )$")
    ax.set_ylabel(r"$n V_{\parallel} / Z \theta_{\parallel b}$")
    ax.set_title(r"$T_{\perp b}/T_{\parallel b}$ = " + str(T_anisotropy))
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.3)
    ax.legend(loc="best", fontsize=15, ncol=2)

plt.tight_layout()
plt.show()