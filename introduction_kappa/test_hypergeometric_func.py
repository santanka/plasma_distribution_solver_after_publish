import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern Roman']
mpl.rcParams['mathtext.fontset'] = 'cm'

plt.rcParams["font.size"] = 20

mp.mp.dps = 50  # 精度（必要なら上げる）

def a_kappa(kappa):
    k = mp.mpf(kappa)
    return (k - 1) / (k - mp.mpf('0.5'))

def I_kappa_hyp2f1(z, kappa):
    """
    I_kappa(z) = (2k/(4k+1)) * 2F1(k+1/2, 1; 2k+3/2; 1 - a_k z)
    """
    k = mp.mpf(kappa)
    ak = a_kappa(k)
    arg = 1 - ak * mp.mpf(z)
    pref = (2*k) / (4*k + 1)
    return pref * mp.hyp2f1(k + mp.mpf('0.5'), 1, 2*k + mp.mpf('1.5'), arg)

def I_limit(z):
    # kappa -> infty の極限（前の議論）
    z = mp.mpf(z)
    return 1 / (1 + z)

def mp_to_float(x):
    # 実数のときだけ float 化（微小虚部は落とす）
    return float(mp.re(x))

# ---- plot settings ----
kappas = [2.0, 5.0, 10.0, 30.0]   # 例。必要に応じて追加
zmin, zmax, Nz = 0.0, 10.0, 400         # arg=1-a_k z が 1 から負側へ動く領域
zs = np.linspace(zmin, zmax, Nz)

fig = plt.figure(figsize=(10, 5), dpi=100)
ax = fig.add_subplot(1, 1, 1, title=r"$I_{\kappa}(z) = \frac{2 \kappa}{4 \kappa + 1} \, {}_{2}F_{1} \left( \kappa + \frac{1}{2}, 1; 2 \kappa + \frac{3}{2}; 1 - \frac{2 (\kappa - 1)}{2 \kappa - 1} z \right)$")

for kappa in kappas:
    ys = []
    for z in zs:
        y = I_kappa_hyp2f1(z, kappa)
        ys.append(mp_to_float(y))
    ax.plot(zs, ys, label=rf"$\kappa$ = {kappa}")

# 極限も参考に
ys_lim = [mp_to_float(I_limit(z)) for z in zs]
ax.plot(zs, ys_lim, linestyle="--", label=r"limit $1/(1+z)$")

ax.set_xlabel(r"$z$")
ax.set_ylabel(r"$I_{\kappa}(z)$")
ax.minorticks_on()
ax.grid(True, alpha=0.3, which="both")
ax.legend()
plt.tight_layout()
plt.show()
