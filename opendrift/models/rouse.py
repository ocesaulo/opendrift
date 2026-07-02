"""Near-bed Rouse resuspension-height sampling for SedimentDrift (Update 2).

Vendored pure-NumPy helpers (no external deps), same self-contained pattern as
``bblm_sg2000.py`` / ``floc_strength.py``. Used by the 'turbulent' resuspension
height mode, which parameterizes the (typically unresolved) near-bed turbulence
directly from the friction velocity ``u*``:

* eddy diffusivity from the law of the wall ``K_z = kappa u* z (1 - z/h)``;
* the equilibrium suspended-sediment distribution is then Rouse. Near the bed
  (z << h) the concentration pdf reduces to a power law ``C(z) ~ z**(-P)`` with
  the Rouse number ``P = w_s / (beta kappa u*)``, whose inverse CDF is
  closed-form -> one cheap, reproducible stochastic draw per particle.

The draw uses the caller-supplied uniform sample (the model passes ``np.random``
so the global, seeded RNG governs reproducibility).
"""

import numpy as np

KAPPA = 0.4  # von Karman constant


def rouse_number(w_s, u_star, beta=1.0, kappa=KAPPA):
    """Rouse number P = w_s / (beta kappa u*). Larger P -> harder to suspend."""
    w_s = np.asarray(w_s, dtype=float)
    u_star = np.asarray(u_star, dtype=float)
    return w_s / np.maximum(beta * kappa * u_star, 1e-12)


def sample_height_powerlaw(P, za, zb, U):
    """Inverse-CDF sample of the near-bed Rouse pdf ``f(z) ~ z**(-P)`` on [za, zb].

    Closed form: for P != 1, ``z = (za**(1-P) + U (zb**(1-P) - za**(1-P)))**(1/(1-P))``;
    for P == 1, ``z = za (zb/za)**U``. All args broadcastable; ``U`` ~ Uniform(0,1).
    Returns the height above the bed, clipped to [za, zb].
    """
    P = np.asarray(P, dtype=float)
    za = np.asarray(za, dtype=float)
    zb = np.maximum(np.asarray(zb, dtype=float), za * (1.0 + 1e-6))
    U = np.asarray(U, dtype=float)

    one_mP = 1.0 - P
    safe = np.where(np.abs(one_mP) < 1e-9, 1.0, one_mP)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        za_p = za ** one_mP
        zb_p = zb ** one_mP
        z_gen = (za_p + U * (zb_p - za_p)) ** (1.0 / safe)
    z_log = za * (zb / za) ** U                      # P -> 1 limit
    z = np.where(np.abs(one_mP) < 1e-6, z_log, z_gen)
    return np.clip(z, za, zb)


def resuspension_height(u_star, w_s, h, za, seed_layer,
                        beta=1.0, P_crit=2.5, U=None):
    """Seed height above the bed [m] for the 'turbulent' mode.

    Particles with Rouse number ``P >= P_crit`` cannot be lofted into suspension
    and are placed in the near-bed (bedload) layer at ``za``. Suspendable
    particles draw a height from the near-bed Rouse power-law profile between
    ``za`` and ``zb = min(h, za + seed_layer)``. ``U`` is a uniform(0,1) sample
    (e.g. ``np.random.uniform`` for the seeded global RNG); if None, drawn here.
    """
    u_star = np.asarray(u_star, dtype=float)
    w_s = np.asarray(w_s, dtype=float)
    h = np.asarray(h, dtype=float)
    za = np.minimum(np.asarray(za, dtype=float), 0.1 * h)   # keep za < h
    # seed_layer may be a per-element array (e.g. the adaptive kappa*u*dt cap).
    zb = np.minimum(h, za + np.asarray(seed_layer, dtype=float))

    n = np.broadcast(u_star, w_s, h, za).shape
    if U is None:
        U = np.random.uniform(size=n)
    P = rouse_number(w_s, u_star, beta=beta)
    z = sample_height_powerlaw(P, za, zb, U)
    return np.where(P >= P_crit, za, z)
