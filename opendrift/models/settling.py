"""Terminal (settling) velocity closures for SedimentDrift.

Vendored alongside ``bblm_sg2000.py`` / ``floc_strength.py`` / ``rouse.py`` so the
installed SedimentDrift model stays self-contained: it does NOT import the offline
research code (``settling_models.py`` / ``maggi.py``, which live outside this
repository in the ``ddt_dump`` project). The relevant
physics is small and reproduced here as pure-NumPy, fully vectorized functions,
with the bugs of those research modules fixed:

* ``settling_models.BB16_shape_factors`` referenced an undefined ``d_eq`` and used
  ``K_S = 0.5*(F_S**(1/3) - F_S**(-1/3))`` (a minus, so a sphere wrongly gave
  ``K_S = 0``). Here the volume-equivalent diameter is handled explicitly and the
  correct ``K_S = 0.5*(F_S**(1/3) + F_S**(-1/3))`` is used (sphere -> 1).
* ``maggi.D0_from_Df`` was a no-op (``D*(D/D)**(Df-3) == D``) and ``maggi_H_approx``
  discarded its result; the fractal-aggregate density used here is the validated
  ``rho_f + (rho_s-rho_f)*(d/d0)**(Df-3)`` (same mass-size relation as
  ``floc_strength.solids_fraction``).

Conventions
-----------
* ``nu`` is the KINEMATIC viscosity [m^2 s^-1] (= dynamic / rho_f). The caller
  (the model) converts the configured molecular dynamic viscosity with the local
  fluid density.
* All wrappers return a POSITIVE settling SPEED [m s^-1]; the caller applies the
  sign (terminal_velocity is negative downward in OpenDrift).
* Every function accepts scalars or NumPy arrays and is vectorized, so one call
  covers all seeded elements.

References: Schiller & Naumann (1933); Dietrich (1982); Bagheri & Bonadonna
(2016); Maggi (2013); Carman-Kozeny; Brinkman.
"""

import numpy as np

GRAV = 9.81  # m s^-2
_TINY = 1e-12


# =========================================================================== #
# Drag laws  (Re = w d / nu, clipped > 0)
# =========================================================================== #
def cd_schiller_naumann(Re):
    """Schiller-Naumann drag coefficient (spheres, Re < ~800)."""
    Re = np.maximum(np.asarray(Re, dtype=float), _TINY)
    return (24.0 / Re) * (1.0 + 0.15 * Re ** 0.687)


def cd_dietrich(Re):
    """Dietrich (1982) drag coefficient (adds the Newton-regime plateau)."""
    Re = np.maximum(np.asarray(Re, dtype=float), _TINY)
    return (24.0 / Re) * (1.0 + 0.15 * Re ** 0.687) + 0.42 / (1.0 + 42500.0 / Re ** 1.16)


def cd_bb16(Re, Ks, Kn):
    """Bagheri & Bonadonna (2016) drag coefficient for non-spherical particles.

    ``Ks``/``Kn`` are the Stokes'/Newton's drag correction factors (1 for a
    sphere); see :func:`bb16_drag_factors`.
    """
    Re = np.maximum(np.asarray(Re, dtype=float), _TINY)
    ReKn = Re * Kn / Ks
    return (24.0 * Ks / Re) * (1.0 + 0.125 * ReKn ** (2.0 / 3.0)) \
        + (0.46 * Kn) / (1.0 + 5330.0 / np.maximum(ReKn, _TINY))


# =========================================================================== #
# BB16 shape -> drag correction factors
# =========================================================================== #
def bb16_drag_factors(flatness, elongation, rho_prime, d_eq=None, L=None, I=None, S=None):
    """Stokes' (Ks) and Newton's (Kn) drag correction factors, Bagheri & Bonadonna (2016).

    ``flatness f = S/I``, ``elongation e = I/L`` from the particle form dimensions
    L >= I >= S. ``rho_prime`` is the particle/fluid density ratio. The form
    factors are

        F_S = f * e**1.3 * (d_eq**3 / (L*I*S))
        F_N = f**2 * e   * (d_eq**3 / (L*I*S))

    where ``d_eq`` is the volume-equivalent sphere diameter. When the true volume
    is unknown (the usual case when only a shape scalar is prescribed) the
    bounding-box estimate ``d_eq**3 = L*I*S`` is used, i.e. the ratio is 1 and
    ``F_S = f*e**1.3``, ``F_N = f**2*e`` (the BB16 form-dimension-only estimate).
    """
    f = np.asarray(flatness, dtype=float)
    e = np.asarray(elongation, dtype=float)
    rho_prime = np.maximum(np.asarray(rho_prime, dtype=float), 1.0 + 1e-9)

    if d_eq is not None and None not in (L, I, S):
        ratio = np.asarray(d_eq, dtype=float) ** 3 / np.maximum(
            np.asarray(L, dtype=float) * np.asarray(I, dtype=float) * np.asarray(S, dtype=float), _TINY)
    else:
        ratio = 1.0

    F_S = np.maximum(f * e ** 1.3 * ratio, _TINY)
    F_N = np.maximum(f ** 2 * e * ratio, _TINY)

    alpha2 = 0.45 + 10.0 / (np.exp(2.5 * np.log(rho_prime)) + 30.0)
    beta2 = 1.0 - 37.0 / (np.exp(3.0 * np.log(rho_prime)) + 100.0)

    Ks = 0.5 * (F_S ** (1.0 / 3.0) + F_S ** (-1.0 / 3.0))          # sphere -> 1
    Kn = 10.0 ** (alpha2 * (-np.log10(F_N)) ** beta2)              # sphere -> 1
    return Ks, Kn


def bb16_factors_from_csf(corey_shape_factor, rho_prime, elongation=1.0):
    """Convenience: BB16 Ks, Kn from a single Corey shape factor.

    The Corey shape factor ``csf = S / sqrt(L*I)`` relates to flatness/elongation
    by ``csf = f * sqrt(e)`` (with ``f = S/I``, ``e = I/L``), so for a prescribed
    elongation (default 1, an equant grain) ``f = csf / sqrt(e)``. csf = 1 (sphere)
    gives Ks = Kn = 1 and BB16 reduces to the sphere drag.
    """
    csf = np.clip(np.asarray(corey_shape_factor, dtype=float), 1e-3, 1.0)
    e = np.asarray(elongation, dtype=float)
    f = np.clip(csf / np.sqrt(e), 1e-3, 1.0)
    return bb16_drag_factors(f, e, rho_prime)


# =========================================================================== #
# Fractal-aggregate (Maggi) density / permeability
# =========================================================================== #
def maggi_effective_density(d, d0, rho_s, rho_f, Df):
    """Effective density of a fractal aggregate (Maggi 2013 mass-size relation).

    ``rho_eff = rho_f + (rho_s - rho_f) * (d/d0)**(Df-3)`` with primary-particle
    size ``d0``, fractal dimension ``Df`` (< 3, so a larger/looser aggregate is
    less dense). Consistent with ``floc_strength.solids_fraction``.
    """
    d = np.asarray(d, dtype=float)
    d0 = np.asarray(d0, dtype=float)
    ratio = (np.maximum(d, _TINY) / np.maximum(d0, _TINY)) ** (np.asarray(Df, dtype=float) - 3.0)
    ratio = np.clip(ratio, 0.0, 1.0)               # aggregate never denser than solid
    return rho_f + (rho_s - rho_f) * ratio


def _carman_kozeny_permeability(d0, eps):
    eps = np.clip(np.asarray(eps, dtype=float), 1e-6, 1.0 - 1e-6)
    return (np.asarray(d0, dtype=float) ** 2 * eps ** 3) / (180.0 * (1.0 - eps) ** 2)


def _brinkman_omega(d, kappa):
    """Brinkman permeability drag-reduction factor Omega in (0, 1]."""
    beta = np.asarray(d, dtype=float) / (2.0 * np.sqrt(np.maximum(kappa, _TINY)))
    beta = np.maximum(beta, 1.0 + 1e-6)            # keep (1 - 1/beta) > 0
    num = 2.0 * beta ** 2 * (1.0 - 1.0 / beta)
    den = 2.0 * beta ** 2 + 3.0 * (1.0 - 1.0 / beta)
    return np.clip(num / den, _TINY, 1.0)


# =========================================================================== #
# Generic iterative terminal-velocity solver (vectorized fixed point)
# =========================================================================== #
def terminal_velocity(d, rho_p, rho_f, nu, cd_func, cd_args=(), drag_scale=1.0,
                      n_iter=15, tol=1e-6):
    """Terminal settling speed from a force balance, by fixed-point iteration.

    ``w = sqrt(4 g d (rho_p - rho_f) / (3 Cd rho_f))`` with ``Cd = drag_scale *
    cd_func(Re, *cd_args)`` and ``Re = w d / nu``, seeded from the Stokes velocity.
    A small fixed iteration cap keeps the cost predictable; sediment converges in
    ~5-10 iterations. Buoyant (rho_p <= rho_f) particles return 0.
    """
    d = np.asarray(d, dtype=float)
    rho_p = np.asarray(rho_p, dtype=float)
    rho_f = np.asarray(rho_f, dtype=float)
    nu = np.asarray(nu, dtype=float)
    drho = rho_p - rho_f

    w = GRAV * drho * d ** 2 / (18.0 * np.maximum(nu * rho_f, _TINY))   # Stokes guess
    w = np.maximum(w, _TINY)
    A = 4.0 * GRAV * d * drho / (3.0 * np.maximum(rho_f, _TINY))

    for _ in range(n_iter):
        Re = np.maximum(w * d / np.maximum(nu, _TINY), _TINY)
        Cd = drag_scale * cd_func(Re, *cd_args)
        w_new = np.sqrt(np.maximum(A / np.maximum(Cd, _TINY), 0.0))
        if np.all(np.abs(w_new - w) <= tol * np.maximum(w, _TINY)):
            w = w_new
            break
        w = w_new

    return np.where(drho > 0.0, w, 0.0)


# =========================================================================== #
# Model wrappers  (positive settling speed [m s^-1])
# =========================================================================== #
def settling_stokes(d, rho_s, rho_f, nu, **_):
    """Stokes settling speed (no drag correction)."""
    d = np.asarray(d, dtype=float)
    drho = np.asarray(rho_s, dtype=float) - np.asarray(rho_f, dtype=float)
    w = GRAV * drho * d ** 2 / (18.0 * np.maximum(np.asarray(nu, dtype=float) * rho_f, _TINY))
    return np.maximum(w, 0.0)


def settling_dietrich(d, rho_s, rho_f, nu, **_):
    """Dietrich (1982) settling speed for a (sub)spherical grain."""
    return terminal_velocity(d, rho_s, rho_f, nu, cd_dietrich)


def settling_bb16(d, rho_s, rho_f, nu, corey_shape_factor=1.0, elongation=1.0, **_):
    """Bagheri & Bonadonna (2016) settling speed for a non-spherical grain."""
    rho_prime = np.asarray(rho_s, dtype=float) / np.maximum(np.asarray(rho_f, dtype=float), _TINY)
    Ks, Kn = bb16_factors_from_csf(corey_shape_factor, rho_prime, elongation=elongation)
    return terminal_velocity(d, rho_s, rho_f, nu, cd_bb16, cd_args=(Ks, Kn))


def settling_maggi(d, rho_s, rho_f, nu, d0=None, fractal_dim=2.0, **_):
    """Maggi (2013) fractal-aggregate settling speed (impermeable)."""
    if d0 is None:
        d0 = d
    rho_eff = maggi_effective_density(d, d0, rho_s, rho_f, fractal_dim)
    return terminal_velocity(d, rho_eff, rho_f, nu, cd_schiller_naumann)


def settling_maggi_permeable(d, rho_s, rho_f, nu, d0=None, fractal_dim=2.0, **_):
    """Maggi (2013) fractal-aggregate settling speed with Brinkman permeability."""
    if d0 is None:
        d0 = d
    rho_eff = maggi_effective_density(d, d0, rho_s, rho_f, fractal_dim)
    eps = 1.0 - np.clip((np.asarray(d, dtype=float) / np.maximum(np.asarray(d0, dtype=float), _TINY))
                        ** (np.asarray(fractal_dim, dtype=float) - 3.0), 0.0, 1.0)
    kappa = _carman_kozeny_permeability(d0, eps)
    omega = _brinkman_omega(d, kappa)
    return terminal_velocity(d, rho_eff, rho_f, nu, cd_schiller_naumann, drag_scale=omega)


# =========================================================================== #
# Master switch
# =========================================================================== #
_MODELS = {
    'stokes': settling_stokes,
    'dietrich': settling_dietrich,
    'bb16': settling_bb16,
    'maggi': settling_maggi,
    'maggi_permeable': settling_maggi_permeable,
}


def settling_velocity(model, d, rho_s, rho_f, nu, **kwargs):
    """Positive settling speed [m s^-1] for ``model``.

    Parameters
    ----------
    model : {'stokes','dietrich','bb16','maggi','maggi_permeable'}
    d : aggregate / grain diameter [m]
    rho_s : solid (or constituent) density [kg m^-3]
    rho_f : fluid density [kg m^-3]
    nu : kinematic viscosity [m^2 s^-1]
    **kwargs : model-specific (``corey_shape_factor``, ``elongation`` for bb16;
        ``d0``, ``fractal_dim`` for the maggi models). Extra keys are ignored, so
        the same per-element property dict can be passed to any model.
    """
    try:
        fn = _MODELS[model]
    except KeyError:
        raise ValueError("Unknown settling model %r; choose from %s"
                         % (model, sorted(_MODELS)))
    return np.maximum(np.asarray(fn(d, rho_s, rho_f, nu, **kwargs), dtype=float), 0.0)
