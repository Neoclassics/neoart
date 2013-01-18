import numpy as np

from geometry import GeometryStrategy
from species import CompositeSpecies
from transport import TransportCoefficients
from arguments import Arguments, DefaultNumericArgs
import core


class Neoart(object):
    """
    Main Neoart object.

    Example
    -------

    Two species (electron and deuterium) with equal temperature and density.
    The ion temperature gradient is non-zero, while the pressure gradient for
    both species is zero.  Fluxes are evaluated in the banana-plateau regime.

    >>> from neoart import *
    >>> e = ('electron', {'density': 5, 'density gradient': 0,
        ...               'temperature': 10, 'temperature gradient': 0,
        ...               'eps':1e-5})
    >>> d = ('deuterium', {'density': 5, 'density gradient': -1,
    ...                   'temperature': 10, 'temperature gradient': 1})
    ...                   'eps':1e-5})
    >>> species = [e, d]

    Use circular geometry
    >>> eps = 1e-4 #  inverse aspect ratio
    >>> geometry = CircularGeometry(minor_radius=eps*1.65, major_radius=1.65,
    ...                             safety_factor=2, magnetic_field=2.5)

    We use the default numerical parameters, but we change the contribution
    >>> numerical = DefaultNumericArgs()
    >>> numerical.set_contribution('banana-plateau')

    >>> n = Neoart(species, geometry, numerical)
    >>> f = n.fluxes()

    Compute the fluxes with parallel electric field
    >>> f = n.flxues(eparr=0.1)
    """

    def __init__(self, species=None, geometry=None, numerical=None, **kwargs):
        self.set_numerical(numerical)
        self.set_geometry(geometry)
        self.set_species(species)
        self.extra_args = kwargs

    def fluxes(self, **kwargs):
        self._geometry.pre(self)
        fluxes = []
        for a in self.get_args(kwargs):
            self._geometry.at_each_rho(a)
            fluxes.append(self._call_neoart(a))
        self._geometry.post(self)

        return np.array(fluxes)

    def set_geometry(self, geometry=None):
        if geometry is None:
            geometry = GeometryStrategy()

        self._geometry = geometry

    def set_numerical(self, numerical=None):
        if numerical is None:
            numerical = DefaultNumericArgs()
        self._numerical = numerical

    def set_species(self, species):
        self._species = CompositeSpecies(species)

    def get_args(self, kwargs={}):
        args = {}
        args.update(self.extra_args)
        args.update(kwargs)
        return Arguments(self._species, self._geometry, self._numerical,
                         extra_args=args)

    def _call_neoart(self, args):
        # Neoart accepts any flux surface label as 'rho', but we prefer to use
        # the more special 'eps' (inverse aspect ratio).  Before calling
        # neoart we rename the key 'eps' to 'rho'.
        args['rho'] = args.pop('eps')
        return core.neoart(**args)

    def get_transport_coeffs(self):
        return TransportCoefficients(self).compute()

    species = property(lambda self: self._species)
    geometry = property(lambda self: self._geometry)
    numerical = property(lambda self: self._numerical)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
