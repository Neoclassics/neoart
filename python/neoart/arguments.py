from collections import Mapping
import numpy as np

from species import CompositeSpecies
from forces import Forces

class Arguments(object):

    def __init__(self, species, geometry, numerical, extra_args=None):
        self.builders = {}
        self.args = {}

        if extra_args is None:
            self._extra_args = {}
        else:
            self._extra_args = extra_args

        self.builders['species'] = SpeciesArgsBuilder(species)
        self.builders['geometry'] = GeometryArgsBuilder(geometry)
        self.builders['numerical'] = NumericArgsBuilder(numerical)

        self.build()

    def build(self):
        for b in self.builders.values():
            b.build()
            self.args.update(b.get_result())

    def get_args(self):
        ret = {}
        ret.update(self.args)
        ret.update(self._extra_args)
        return ret

    def __iter__(self):
        iterables, args = self._split_iterables()
        for v in zip(*iterables.values()):
            kw = dict(zip(iterables.keys(), v))
            kw.update(args)
            kw.update(self._extra_args)
            yield kw

    def _split_iterables(self):
        args = self.get_args().copy()
        keys = ['forces', 'temperature', 'density', 'eps']

        iterables = {}
        for k in keys:
            iterables[k] = args.pop(k)

        return iterables, args


# Builder interface for the input arguments
class ArgsBuilder(object):

    argnames = ''

    def __init__(self):
        self.args = {}

    def build(self):
        pass

    def get_result(self):
        return self.args

    def delivers(self):
        return [t.strip() for t in self.argnames.split(',')]


class GeometryArgsBuilder(ArgsBuilder):

    argnames = 'isel'

    def __init__(self, geometry):
        super(GeometryArgsBuilder, self).__init__()
        self._geometry = geometry

    def build(self):
        self.args.update(self._geometry.get_args())


class NumericArgsBuilder(ArgsBuilder):

    argnames = 'ic, force_viscosity, electron_ion_collisions'

    def __init__(self, data=None):
        super(NumericArgsBuilder, self).__init__()
        self.contribution = None
        self.eparr = 0

        if isinstance(data, NumericArgsBuilder):
            self.set_contribution(data.contribution)
            self.set_electron_ion_collisions(data.electron_ion_collisons)
            self.set_force_viscosity_regime(data.force_viscosity_regime)
            self.set_eparr(data.eparr)
        elif isinstance(data, Mapping):
            for key, value in data.iteritems():
                attr = getattr(self, 'set_%s' % key)
                attr(value)

    def set_contribution(self, value):
        if value in ['cl', 'classical']:
            self.contribution = 'classical'
        elif value in ['bp', 'banana-plateau']:
            self.contribution = 'banana-plateau'
        elif value in ['ps', 'pfirsch-schlueter']:
            self.contribution = 'pfirsch-schlueter'
        elif value in ['all']:
            self.contribution = 'all'
        else:
            raise NotImplementedError('contribution %s is unknown.'
                                      % value)

    def set_electron_ion_collisions(self, value):
        self.electron_ion_collisons = value

    def set_force_viscosity_regime(self, value):
        self.force_viscosity_regime = value

    def set_eparr(self, value):
        self.eparr = value

    def _build_contribution(self):
        ic_codes = {'classical': 0, 'banana-plateau': 1,
                    'pfirsch-schlueter': 2, 'all': 3}
        return ic_codes[self.contribution]

    def build(self):
        self.args.update({'ic': self._build_contribution()})
        self.args.update({'force_viscosity': self.force_viscosity_regime})
        self.args.update({'electron_ion_collisions':
                          self.electron_ion_collisons})
        self.args.update({'eparr': self.eparr})


class DefaultNumericArgs(NumericArgsBuilder):
    """
    Default numerical arguments.
    """

    def __init__(self, data=None):
        super(DefaultNumericArgs, self).__init__(data)
        self.set_contribution('all')
        self.set_electron_ion_collisions(True)
        self.set_force_viscosity_regime(0)
        self.set_eparr(0)


class SpeciesArgsBuilder(ArgsBuilder):

    argnames = 'eps, nc, mass, zsp, density, temperature, forces'

    def __init__(self, species):
        super(SpeciesArgsBuilder, self).__init__()
        self._species = CompositeSpecies(species)

    def build(self):
        zsp = self._species.get_zsp('density')
        nc = [len(x) for x in zsp]

        self.args.update({'mass': self._species.get_mass()})
        self.args.update({'nc': nc})
        self.args.update({'zsp': np.array(np.broadcast_arrays(*zsp))})
        self.args.update(self._build_density_temperature())
        self.args.update(self._build_forces())

    def _build_density_temperature(self):
        d = {}
        for p in ['density', 'temperature', 'eps']:
            profiles = self._species.get_profile(p)
            profiles = np.array(np.broadcast_arrays(*profiles))
            d[p] = myrollaxis(profiles, 3)

        d['temperature'] = self.remove_dim_charge_state(d['temperature'])
        d['eps'] = self.remove_dim_charge_state(d['eps'])[..., 0]

        return d

    def _build_forces(self):
        gradients = []
        for s in self._species:
            forces = self.get_forces(s)
            pressure_gradient = forces.get_pressure_gradient()
            temperature_gradient = forces.get_temperature_gradient()
            gradients.append(
                np.concatenate((pressure_gradient[..., None, :],
                                temperature_gradient[..., None, :]), axis=-2))

        forces = np.array(np.broadcast_arrays(*gradients))
        return {'forces': myrollaxis(forces, 4)}

    def get_forces(self, s):
        return Forces(s.profiles, squeeze=False)

    def remove_dim_charge_state(self, a):
        if a.ndim == 3:
            return a[..., 0]
        else:
            return a


def myrollaxis(a, min_ndim):
    if a.ndim == min_ndim:
        return np.rollaxis(a, min_ndim - 1)
    elif a.ndim == min_ndim - 1:
        return a[None]


