"""Critical-shear-stress closures for SedimentDrift resuspension.

Vendored alongside ``bblm_sg2000.py`` so the installed SedimentDrift model stays
self-contained: it does NOT import the offline ``codes/maggi.py`` /
``codes/settling_models.py``. The relevant physics is small and reproduced here
as pure-NumPy, fully vectorized functions.

Two regimes (see RESUSPENSION_MECHANICS_PLAN.md, Update 1):

* Non-cohesive grains -- Soulsby & Whitehouse (1997) fit to the Shields curve,
  a function of the grain's own size/density and the ambient fluid.
* Cohesive flocs -- Kranenburg (1994) / Winterwerp & van Kesteren fractal bed
  strength: the erosion threshold scales as ``phi**(2/(3-Df))`` (phi = solids
  volume fraction of the aggregate, from the Maggi 2013 mass-size relation)
  times a Shields-like buoyant-stress scale ``g*(rho_s-rho_f)*d_floc``.
"""

import numpy as np

GRAV = 9.81  # m s^-2


def solids_fraction(d_floc, d0, Df):
    """Solids volume fraction ``phi`` of a fractal aggregate.

    From the Maggi (2013) mass-size relation ``phi = (d_floc/d0)**(Df-3)`` with
    ``d0`` the primary-particle size and ``Df`` the fractal dimension (Df<3, so a
    larger/looser floc has a smaller solids fraction). Clipped to (0, 1].
    """
    d_floc = np.asarray(d_floc, dtype=float)
    d0 = np.asarray(d0, dtype=float)
    Df = np.asarray(Df, dtype=float)
    phi = (np.maximum(d_floc, 1e-12) / np.maximum(d0, 1e-12)) ** (Df - 3.0)
    return np.clip(phi, 1e-6, 1.0)


def consolidate_phi(phi0, t_settled, t_consol, phi_max):
    """Solids fraction after consolidation.

    Relaxes from the fresh value ``phi0`` toward a packed ``phi_max`` over the
    consolidation timescale ``t_consol`` (Sanford 2008-style dewatering):
    ``phi(t) = phi0 + (phi_max - phi0) * (1 - exp(-t_settled / t_consol))``.
    """
    phi0 = np.asarray(phi0, dtype=float)
    phi_max = np.asarray(phi_max, dtype=float)
    t_settled = np.asarray(t_settled, dtype=float)
    f = 1.0 - np.exp(-t_settled / max(float(t_consol), 1e-9))
    phi = phi0 + (np.maximum(phi_max, phi0) - phi0) * f
    return np.clip(phi, 1e-6, 1.0)


def tau_crit_floc(phi, Df, d_floc, rho_s, rho_f, c_str=1.0):
    """Cohesive critical shear stress [Pa] from fractal floc strength.

    ``tau = c_str * g*(rho_s-rho_f)*d_floc * phi**(2/(3-Df))`` -- a Shields-like
    buoyant-stress scale modulated by the fractal cohesion factor. Depends only
    on the particle's own floc state (d_floc, rho_s, Df, phi) and the fluid rho_f.
    """
    phi = np.clip(np.asarray(phi, dtype=float), 1e-6, 1.0)
    Df = np.asarray(Df, dtype=float)
    d_floc = np.asarray(d_floc, dtype=float)
    rho_s = np.asarray(rho_s, dtype=float)
    rho_f = np.asarray(rho_f, dtype=float)
    expo = 2.0 / np.maximum(3.0 - Df, 1e-3)
    drho = np.maximum(rho_s - rho_f, 0.0)
    return c_str * GRAV * drho * d_floc * phi ** expo


def tau_crit_shields(d, rho_s, rho_f, nu_dyn):
    """Non-cohesive critical shear stress [Pa] (Soulsby & Whitehouse 1997 fit).

    ``theta_cr = 0.30/(1+1.2 D*) + 0.055 (1 - exp(-0.020 D*))`` with the
    dimensionless grain size ``D* = d (g(s-1)/nu^2)^(1/3)``, ``s = rho_s/rho_f``;
    ``tau = theta_cr * g (rho_s-rho_f) d``. ``nu_dyn`` is the dynamic (molecular)
    viscosity [kg m^-1 s^-1]; converted to kinematic with the local rho_f.
    """
    d = np.asarray(d, dtype=float)
    rho_s = np.asarray(rho_s, dtype=float)
    rho_f = np.asarray(rho_f, dtype=float)
    nu_kin = np.asarray(nu_dyn, dtype=float) / rho_f          # m^2 s^-1
    s = rho_s / rho_f
    Dstar = d * (GRAV * np.maximum(s - 1.0, 1e-6) / nu_kin ** 2) ** (1.0 / 3.0)
    theta_cr = 0.30 / (1.0 + 1.2 * Dstar) + 0.055 * (1.0 - np.exp(-0.020 * Dstar))
    return theta_cr * GRAV * (rho_s - rho_f) * d
