"""
Styles & Glenn (2000) neutral combined wave-current bottom boundary layer model.

Vectorized NumPy port of the reference MATLAB implementation in this folder
(bblm02.m / bstress2.m / phi2_1.m / pwave.m / shldc.m). All quantities are in
**cgs units** internally (cm, cm/s, cm/s^2), matching the MATLAB code and the
UserTestCase. The integration layer (SedimentDrift) is responsible for any
SI<->cgs conversion.

Inputs are the near-bed wave + current state, supplied generically so they can be
sourced from a wave hindcast reader (ERA5/WW3/CDIP), a wind parameterisation, or a
test fixture:
    Ub  - near-bed wave orbital velocity      [cm/s]
    Ab  - near-bed wave excursion amplitude    [cm]   (wave radian freq omega = Ub/Ab)
    Ur  - mean current speed at height zr       [cm/s]
    zr  - height above bed where Ur is given     [cm]
    deg - angle between wave and current          [degrees]

Reference:
    Styles, R. and S. M. Glenn, 2000: Modeling stratified wave and current bottom
    boundary layers on the continental shelf. JGR, 105(C10), 24119-24139.
    Styles, Glenn & Brown, 2017: ERDC/CHL TR-17-11 (optimized arbitrary roughness).

Validation: see test_sg2000.py (reproduces UserTestCase Table 2 from test2.mat).
"""

import numpy as np
from scipy.special import kelvin


# ----------------------------------------------------------------------
# Default physical / closure constants (cgs), from bblm02.m
# ----------------------------------------------------------------------
DEFAULTS = dict(
    d_median=0.04,   # median grain diameter [cm]
    s=2.65,          # relative sediment density
    nu=0.0119,       # kinematic viscosity of seawater @15C [cm^2/s]
    g=981.0,         # gravity [cm/s^2]
    kappa=0.4,       # von Karman constant
    beta=0.7,        # closure constant
    Alpha=0.3,       # closure constant (0.3 recommended, sandy bed w/ waves)
    Con=6.4,         # ripple-roughness coefficient paired with Alpha=0.3
    eta_def=1.0,     # default ripple height [cm] (no sediment motion)
    lam_def=15.0,    # default ripple wavelength [cm]
    kbr_def=3.0,     # default ripple roughness [cm]
    tol=1.0e-4,
    max_bisect=50,   # mm in bblm02.m
    max_secant=40,   # pwave.m iteration cap
)


def shields_critical(star):
    """Critical Shields parameter (initiation of motion) from the Shields diagram.

    Vectorized port of shldc.m. `star` is the fluid-sediment parameter
    star = d/(4 nu) * sqrt(g d (s-1)).
    """
    star = np.asarray(star, dtype=float)
    conds = [star < 1.5, star < 4.0, star < 10.0, star < 34.0, star < 270.0]
    choices = [
        0.0932 * star ** (-0.707),
        0.0848 * star ** (-0.473),
        0.0680 * star ** (-0.314),
        np.full_like(star, 0.033),
        0.0134 * star ** 0.255,
    ]
    return np.select(conds, choices, default=np.full_like(star, 0.056))


def _kelvin_phi(znotp, z1p, mp, kappa):
    """Non-dimensional wave shear phi (vectorized port of phi2_1.m / pwave inner loop).

    Two-layer Kelvin-function solution where z1p/znotp > 1, else the analytic
    limit |-kappa z1p mp|. scipy.special.kelvin(x) returns
    (ber+i bei, ker+i kei, ber'+i bei', ker'+i kei').
    """
    znotp = np.asarray(znotp, dtype=float)
    z1p = np.broadcast_to(np.asarray(z1p, dtype=float), znotp.shape)
    mp = np.broadcast_to(np.asarray(mp), znotp.shape)

    # analytic limit (used where z1p/znotp <= 1, and as a safe default)
    phi_analytic = np.abs(-kappa * z1p * mp)

    use_kelvin = (z1p / znotp) > 1.0
    phi = phi_analytic.copy()
    if np.any(use_kelvin):
        zn = znotp[use_kelvin]
        z1 = z1p[use_kelvin]
        mpk = mp[use_kelvin]

        Be_n, Ke_n, Bep_n, Kep_n = kelvin(2.0 * np.sqrt(zn))
        Be_1, Ke_1, Bep_1, Kep_1 = kelvin(2.0 * np.sqrt(z1))

        bnot = Be_n
        knot = Ke_n
        bnotp = Bep_n / np.sqrt(zn)
        knotp = Kep_n / np.sqrt(zn)

        b1 = Be_1
        k1 = Ke_1
        b1p = Bep_1 / np.sqrt(z1)
        k1p = Kep_1 / np.sqrt(z1)

        ll = mpk * b1 + b1p
        nn = mpk * k1 + k1p
        denom = bnot * nn - knot * ll
        argi = bnotp * nn / denom - knotp * ll / denom
        gammai = -kappa * zn * argi
        phi[use_kelvin] = np.abs(gammai)
    return phi


def _residual(uboucw, abozn, zrozn, ubokur, theta, alpha_cl, kappa, z1p, mp):
    """f(uboucw) whose root gives uboucw = Ub/u*cw. Port of bstress2.m.

    Returns (fofx, info) where info holds Ro, mu, epsilon, z1ozn, z2ozn,
    zroz1, zroz2 for the converged solution.
    """
    uboucw = np.asarray(uboucw, dtype=float)
    Ro = abozn / uboucw
    znotp = 1.0 / (kappa * Ro)

    phi = _kelvin_phi(znotp, z1p, mp, kappa)
    mu = np.sqrt(uboucw * phi)

    # epsilon = u*c / u*cw  (clamp eps^2 >= 0 for robustness, cf. commented MATLAB fallback)
    eps2 = -mu ** 2 * np.abs(np.cos(theta)) + np.sqrt(
        np.maximum(1.0 - mu ** 4 * np.abs(np.sin(theta) ** 2), 0.0))
    eps2 = np.maximum(eps2, 1.0e-12)
    epsilon = np.sqrt(eps2)

    Ror = Ro / zrozn
    zroz1 = 1.0 / (alpha_cl * kappa * Ror)
    zroz2 = epsilon * zroz1
    z1ozn = alpha_cl * kappa * Ro
    z2ozn = z1ozn / epsilon

    with np.errstate(invalid="ignore", divide="ignore"):
        f1 = ubokur * epsilon * (np.log(zroz2) + 1 - epsilon + epsilon * np.log(z1ozn)) - uboucw
        f2 = ubokur * epsilon ** 2 * (zroz1 - 1 + np.log(z1ozn)) - uboucw
        f3 = ubokur * epsilon ** 2 * np.log(zrozn) - uboucw
        f4 = ubokur * epsilon * (np.log(zroz2) + 1 - 1.0 / z2ozn) - uboucw
        f5 = ubokur * epsilon ** 2 * (zroz1 - 1.0 / z1ozn) - uboucw
        f6 = ubokur * epsilon * np.log(zrozn) - uboucw  # same form as f3

    conds = [
        (zroz2 > 1) & (z1ozn > 1),
        (zroz2 <= 1) & (zroz1 > 1) & (z1ozn > 1),
        (zroz1 <= 1) & (z1ozn > 1),
        (zroz2 > 1) & (z1ozn <= 1) & (z2ozn > 1),
        (zroz2 <= 1) & (zroz1 > 1) & (z1ozn <= 1) & (z2ozn > 1),
        (zroz2 > 1) & (z2ozn <= 1),
    ]
    fofx = np.select(conds, [f1, f2, f3, f4, f5, f6], default=f1)

    info = dict(Ro=Ro, mu=mu, epsilon=epsilon, z1ozn=z1ozn,
                z2ozn=z2ozn, zroz1=zroz1, zroz2=zroz2)
    return fofx, info


def _pure_wave_limit(abozn, alpha_cl, kappa, z1p, mp, tol, max_secant):
    """Pure-wave (current -> 0) limit of uboucw, an upper bracket. Port of pwave.m.

    Vectorized secant iteration on f(u) = u - 1/phi(u).
    """
    abozn = np.asarray(abozn, dtype=float)
    z1p = np.broadcast_to(np.asarray(z1p, dtype=float), abozn.shape)
    mp = np.broadcast_to(np.asarray(mp), abozn.shape)
    t1 = -alpha_cl * kappa * z1p

    # initial guess (piecewise on abozn), as in bblm02.m
    logb = np.log(np.maximum(abozn, 1e-30))
    g0 = 1.0 / np.abs(t1 * mp)
    g1 = np.exp(1.488) * abozn ** (-0.653) * abozn ** (0.185 * logb)
    g2 = np.exp(0.4599) * abozn ** (0.1977) * abozn ** (0.0085 * logb)
    g3 = np.exp(0.13996) * abozn ** (0.3539) * abozn ** (-0.0106 * logb)
    ubouwmgs = np.select(
        [abozn < 6.25, abozn < 10.0, abozn < 100.0],
        [g0, g1, g2], default=g3)

    # secant (replicates pwave.m bookkeeping)
    ubouwm = ubouwmgs.copy()
    ubouwmo = ubouwmgs * 0.293847
    ubouwmn = 0.5234511947 * ubouwmgs
    phio = np.full_like(abozn, 0.9)

    def phi_of(u):
        u = np.where(u < 0, 1.0e-8, u)
        Ro = abozn / u
        znotp = 1.0 / (kappa * Ro)
        return _kelvin_phi(znotp, z1p, mp, kappa)

    active = np.ones_like(abozn, dtype=bool)
    for cnt in range(max_secant):
        if cnt > 0:
            ubouwmo = np.where(active, ubouwm, ubouwmo)
            ubouwm = np.where(active, ubouwmn, ubouwm)
            phio = np.where(active, phi, phio)
        ubouwm = np.where(ubouwmn < 0, 1.0e-8, ubouwm)
        phi = phi_of(ubouwm)
        fofsigma = ubouwm - 1.0 / phi
        fofsigmao = ubouwmo - 1.0 / phio
        denom = (fofsigma - fofsigmao)
        denom = np.where(denom == 0, 1e-30, denom)
        new = ubouwm - fofsigma * (ubouwm - ubouwmo) / denom
        converged = np.abs((new - ubouwm) / np.where(new == 0, 1e-30, new)) <= tol
        ubouwmn = np.where(active, new, ubouwmn)
        active = active & ~converged
        if not np.any(active):
            break
    return ubouwmn


def sg2000_solve(Ub, Ab, Ur, zr, deg, **kw):
    """Solve the SG2000 neutral wave-current BBL for shear velocities.

    All inputs cgs; scalars or broadcastable arrays. Returns a dict of arrays:
        ustarcw - combined max wave-current shear velocity [cm/s]
        ustarc  - time-mean current shear velocity          [cm/s]
        ustarwm - max wave shear velocity                    [cm/s]
        znot    - hydraulic roughness                         [cm]
        z1, z2  - inner/transition BBL heights                [cm]
        mu, epsilon, Ro, uboucw  - non-dimensional diagnostics
        kbr, kbs - ripple / sediment-transport roughness      [cm]
    """
    p = {**DEFAULTS, **kw}
    d_median, s, nu, g = p["d_median"], p["s"], p["nu"], p["g"]
    kappa, beta, Alpha, Con = p["kappa"], p["beta"], p["Alpha"], p["Con"]

    Ub, Ab, Ur, zr, deg = np.broadcast_arrays(
        *[np.asarray(x, dtype=float) for x in (Ub, Ab, Ur, zr, deg)])
    theta = deg * np.pi / 180.0
    omega = Ub / Ab

    # skin-friction Shields parameter (Madsen) vs critical (Shields diagram)
    psinorm = (s - 1) * g * d_median
    arg_ole = Ab / d_median
    fwcskn = np.exp(5.61 * arg_ole ** (-0.109) - 7.30)
    Psi = 0.5 * fwcskn * (1.42 * Ub) ** 2 / psinorm
    star = d_median / (4 * nu) * np.sqrt(g * d_median * (s - 1))
    psicr = shields_critical(star)

    # ripple geometry from mobility number CHI (SG2002)
    CHI = 4 * nu * Ub ** 2 / (d_median * ((s - 1) * g * d_median) ** 1.5)
    eta = np.where(CHI < 2, Ab * 0.30 * CHI ** (-0.39), Ab * 0.45 * CHI ** (-0.99))

    # bottom roughness: kb = d_median + ripple + sediment-transport roughness
    moves = (Psi - psicr) > 0
    kbr = np.where(moves, Con * eta, p["kbr_def"])
    kbs = Ab * 0.0655 * (Ub ** 2 / ((s - 1) * g * Ab)) ** 1.4
    kb = d_median + kbr + kbs
    znot = kb / 30.0

    abozn = Ab / znot
    zrozn = zr / znot
    ubokur = (Ub / Ur) / kappa

    alpha_cl = Alpha * (1 + beta * kb / Ab)   # closure scale; z1' = alpha_cl
    z1p = alpha_cl
    delta = 1.0 / np.sqrt(2 * z1p)
    mp = delta + 1j * delta

    # bracket [a, b]: a ~ 0, b = pure-wave limit
    a = np.full_like(abozn, 1.0e-6)
    b = _pure_wave_limit(abozn, alpha_cl, kappa, z1p, mp, p["tol"], p["max_secant"])

    fa, _ = _residual(a, abozn, zrozn, ubokur, theta, alpha_cl, kappa, z1p, mp)
    # orient so that f(lo) and f(hi) bracket the root; f decreases with uboucw,
    # so f(a) > 0 and f(b) < 0 in the normal case.
    sign_a = np.sign(fa)

    lo, hi = a.copy(), b.copy()
    c = 0.5 * (lo + hi)
    for _ in range(p["max_bisect"]):
        fc, info = _residual(c, abozn, zrozn, ubokur, theta, alpha_cl, kappa, z1p, mp)
        same_side = np.sign(fc) == sign_a
        lo = np.where(same_side, c, lo)
        hi = np.where(same_side, hi, c)
        c = 0.5 * (lo + hi)
        if np.all((hi - lo) < p["tol"]):
            break
    uboucw = c
    _, info = _residual(uboucw, abozn, zrozn, ubokur, theta, alpha_cl, kappa, z1p, mp)

    Ro = info["Ro"]
    ustarcw = Ro * znot * omega
    ustarc = info["epsilon"] * ustarcw
    ustarwm = info["mu"] * ustarcw

    return dict(
        ustarcw=ustarcw, ustarc=ustarc, ustarwm=ustarwm,
        znot=znot, z1=info["z1ozn"] * znot, z2=info["z2ozn"] * znot,
        mu=info["mu"], epsilon=info["epsilon"], Ro=Ro, uboucw=uboucw,
        kbr=kbr, kbs=kbs, psicr=psicr, Psi=Psi,
    )
