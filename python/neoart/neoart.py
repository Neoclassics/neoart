import numpy as np
import _neoart

from configuration import Configuration
from result import Result


def neoart(rho, config):
    """
    Low level wrapper of NEOART.
    """
    d = config.todict()
    n_species = len(d['masses'])

    shape_ = _neoart.get_elem_config()
    zsp = np.zeros(shape_)
    den = np.zeros(shape_)
    ds = np.zeros(shape_ + (2,))

    den[:n_species, :] = d['den']
    zsp[:n_species, :] = d['zsp']
    ds[:n_species, :, :] = d['ds']

    coeff = _neoart.neoart(d['nc'], zsp, d['masses'], d['temperatures'], den,
            ds, rho, d['eps'], d['isel'], d['ishot'], d['nreg'], d['sigma'],
            d['nleg'], d['nenergy'], d['ncof'], d['neogeo'],
            d['neofrc'], d['ic'], d['eparr'])

    coeff = coeff[:n_species, 0, :]
    return Result(coeff, config)


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


def colxi(config):
    """
    Calculates the collision frequencies and the weighting factors xi.

    Parameters
    ----------
    tau : array(nsm,nsm)
        Collision frequency weighted over charge states.
    xi : array(nsm,nsm)
        Relative weight of every charge state.
    """
    d = config.todict()
    tau, xi = _neoart.colxi(d['nc'], d['zsp'], d['den'], d['temperatures'],
            d['masses'])

    return tau, xi


if __name__ == '__main__':
    import doctest
    doctest.testmod()
