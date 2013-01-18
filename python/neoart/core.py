__all__ = ['neoart', 'colxi']

import numpy as np
import _neoart


def neoart(nc, zsp, mass, temperature, density, forces, rho, isel, eparr=0,
           ic=3, accuracy=1e-5, force_viscosity=0, sigma=(1, 1, 1, 1), nleg=3,
           energy_scattering=True, electron_ion_collisions=True):
    """
    Low-level interface for computing neoclassical fluxes with NEOART.  This
    routine is usually called by the Neoart object which exposes a simpler and
    more intuitive interface.  Usually the user should not need to bother with
    this.

    Parameters
    ----------
    nc : sequence of ints
        Number of charge states per species.  The length of `nc` defines the
        number of species: ns = len(nc)
    zsp : ndarray
        Charge number of all the charge states of all species.
        shape = (number of species, maximal number of charge states)
    mass : array_like
        Species mass in kg. The length of the array must be equal to the
        number of species: len(mass) = nc.
    temperature : array_like
        Species temperature in keV.  The length of the array must be equal to
        the number of species: len(mass) = nc.
    density : array_like
        Density of the charge states in 1e19 m-3.  The shape must be (ns,
        max(nc)).
    forces : array_like
        Pressure and temperature gradients. The shape is (ns, max(nc), 2).

            forces[i, j, 1] = -d ln(p_i_j)/drho
            forces[i, j, 2] = -d ln(T_i)/drho

        for the j-th charge state of the i-th species. ds > 0 for peaked
        profile.
    rho : float
        Flux surface label
    isel : int
        Geometry selection.
    eparr : float
        Loop voltage [V] divided by 2pi.
    ic : int
        The contribution for which the coefficients are calculated.
    accuracy : float
        Accuracy.
    force_viscosity : int
        Force a certain regime in the calculation of the viscosity.
            0     : the default value. the banana and Pfirsch-Schlueter
                    contributions are weighted.
            1     : the banana regime is forced.
            other : the Pfirsch-Schlueter regime is forced.

        Note: forcing the banana regime is not the same thing as calculating
        the banana plateau contribution. for high collision freq. this contri-
        bution usually decreases. when force_viscosity = 1 there is no such
        decr.
    sigma : array_like
        The sigma's determine which of the coupling terms in the equation for
        the Pfirsch-Schlueter regime are taken into account.
    nleg : int
        number of legendre polynomals in the expansion (maximum and normal
        value is 3)
    energy_scattering : bool
        Parameter that determines whether energy scattering is taken into
        account in the calculation of the viscosity.
            True  : accounted for (default value)
            False : not accounted for
    electron_ion_collisions : bool
        Parameters that determines whether ion-electron collisions are taken
        into account or not taken into account (default value: True)

        Note: use False only to obtain some analytic results.

    Returns
    -------
    fluxes : ndarray
        (ns, max(nc), 4) shaped array with the fluxes:
        * particle flux, positive values for outward flux:

            fluxes[i, j, 0] = <Gamma_ij . grad_rho>,

          with Gamma_ij in 1/m^2/s and <.> the flux surface average
        * heat flux, positive values for outward flux

            fluxes[i, j, 1] = <Q_ij . grad_rho>/T_i

          with Q_ij/T_i in 1/m^2/s
        * parallel flow

            fluxes[i, j, 2] = q_ij*n_ij*<u_ij . B><B>/<B^2>

          with q_ij*n_ij*u_ij the current density of species i,j in A/m^2
        * poloidal flow

            fluxes[i, j, 3] = <u_ij.grad_theta>/<B.grad_theta>

          in m/s/T with u_ij the flow of species i,j in m/s.
    """

    nc, m, temp = np.atleast_1d(nc, mass, temperature)
    den, zsp = np.atleast_2d(density, zsp)
    ds = np.atleast_3d(forces)

    n_species = np.alen(nc)
    max_nc = max(nc)

    assert den.shape == (n_species, max_nc)
    assert zsp.shape == (n_species, max_nc)
    assert ds.shape == (n_species, max_nc, 2)

    shape_ = _neoart.get_elem_config()
    zsp_ = np.zeros(shape_)
    den_ = np.zeros(shape_)
    ds_ = np.zeros(shape_ + (2,))

    den_[:n_species, :max_nc] = den
    zsp_[:n_species, :max_nc] = zsp
    ds_[:n_species, :max_nc] = ds

    nreg = force_viscosity
    sigma = np.array(sigma)
    nenergy = energy_scattering
    ncof = electron_ion_collisions

    # Other parameters needed to call the main routine.
    ishot = 0
    neofrc = 0
    neogeo = 1

    coeff = _neoart.neoart(nc, zsp_, mass, temperature, den_, ds_, rho,
                           accuracy, isel, ishot, nreg, sigma, nleg, nenergy,
                           ncof, neogeo, neofrc, ic, eparr)

    coeff[:, :, 0] *= -1.
    coeff[:, :, 1] *= -1.

    coeff = coeff[:n_species, :max_nc, :]

    return np.array(coeff)


def colxi(nc, zsp, mass, temperature, density):
    """
    Calculates the collision frequencies and the weighting factors xi.

    Parameters
    ----------
    tau : array(nsm,nsm)
        Collision frequency weighted over charge states.
    xi : array(nsm,nsm)
        Relative weight of every charge state.
    """
    tau, xi = _neoart.colxi(nc, zsp, density, temperature, mass)
    return tau, xi


def circgeom(rho, Rmag, eps, q, bn):
    """
    Store the parameters of the circular geometry.

    Parmeters
    ---------
    rho : scalar, float
        Flux surface label.
    rn: scalar, float
        Major radius of the magnetic axis.
    e: scalar, float
        Inverse aspect ratio of the surface.
    q: scalar, float
        Safety factor of the flux surface.
    bn: scalar, float
        Magnetic field strentgh in chi = pi / 1 with chi being the poloidal
        angle.

    Example
    -------
    >>> import neoart
    >>> neoart.circgeom(rho=1.65e-4, Rmag=1.65, eps=1e-4, q=2, bn=2.5)
    """
    _neoart.circgeom(1, rho, Rmag, eps, q, bn)
