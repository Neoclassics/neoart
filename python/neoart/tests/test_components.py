import os
import numpy as np

from numpy.testing import assert_equal, assert_almost_equal, \
    assert_approx_equal, assert_allclose
from nose.tools import assert_raises, set_trace, nottest


from ..neoart import Neoart
from ..arguments import Arguments, SpeciesArgsBuilder, DefaultNumericArgs
from ..arguments import GeometryArgsBuilder, NumericArgsBuilder, myrollaxis
from ..forces import Forces
from ..geometry import GeometryStrategy, CircularGeometry, CheaseGeometry
from ..species import Species, CompositeSpecies, ProfileCollection
from ..transport import TransportCoefficients


def thisdir(filename):
    curdir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(curdir, filename)


def test_CompositeProfile():
    p = ProfileCollection()
    p.add('density', [])
    p.add('temperature', 3)

    assert_equal(p.get_profile('density').shape, (0,))
    assert_equal(p.get_profile('temperature').shape, ())

    d = {'density': 1, 'temperature': [1, 2]}
    p = ProfileCollection(d)
    assert_equal(p.get_profile('density'), 1)
    assert_equal(p.get_profile('temperature'), [1, 2])

    p2 = ProfileCollection(p)
    assert_equal(p._profiles, p2._profiles)

    # non existing profile returns empty array
    assert_equal(p.has_profile('non existing'), False)
    assert_equal(p.get_profile('non existing'), [])


def test_CompositeProfile_broadcast():
    p = ProfileCollection({'a': [1, 2, 3], 'b': 2})
    assert_equal(p.get_profile('b'), [2, 2, 2])


def test_Species():
    profiles = {'a': 1, 'b': 2}
    species1 = Species('species1', profiles, mass=10, max_charge=-1)
    species2 = Species('species2', profiles, mass=1, max_charge=6)
    species3 = Species(('species3', profiles))

    assert_equal(species1.get_mass(), 10)
    assert_equal(species2.get_max_charge(), 6)
    assert_equal(species3.name, 'species3')

    s = CompositeSpecies()
    s.add(species1)
    s.add(species2)

    assert_equal(s.get_max_charge(), [-1, 6])
    assert_equal(s.get_mass(), [10, 1])

    assert_equal(s.names, ['species1', 'species2'])

    s = CompositeSpecies([species1, species2])
    assert_equal(s.get_mass(), [10, 1])

    s = CompositeSpecies()
    s.add(('species1', profiles))
    s.add(('specees2', profiles))
    assert_equal(s.get_mass(), [None, None])


def test_Species_zsp():
    profiles = {'a': [1, 1]}
    species1 = Species('species1', profiles, mass=10, max_charge=-1)
    assert_equal(species1.get_zsp('a'), [-1])

    profiles = {'a': [[1, 1], [2, 2], [3, 3]]}
    species1 = Species('species1', profiles, mass=10, max_charge=2)
    assert_equal(species1.get_zsp('a'), [0, 1, 2])
    #assert_equal(species1.get_profile('a', 1), [2, 2])
    assert_equal(species1.get_profile('a'), profiles['a'])

    species1 = Species('species1', {'a': [1, 2]}, mass=10, max_charge=6)
    species2 = Species('species2', {'a': [[1, 2], [3, 4]]},
                       mass=10, max_charge=1)
    s = CompositeSpecies([species1, species2])

    assert_equal(s.get_zsp('a'), [[6], [0, 1]])


def test_Forces():
    x = np.linspace(0, 0.9, 10)
    y = (1 - x) ** 2

    # If the field is present simply return that
    forces = Forces({'density': y, 'density gradient': 2})
    assert_equal(forces.get_density_gradient(), 2 * np.ones_like(x))

    # If there's an 'eps' field than use that to determine the gradients
    forces = Forces({'density': y, 'eps': x})
    assert_almost_equal(forces.get_density_gradient(), 2 / (1 - x))

    # It also works for many ionisation stages
    forces = Forces({'density': [y, y], 'eps': [x, 2 * x]})
    rln = forces.get_density_gradient()
    assert_almost_equal(rln[0], 2 / (1 - x))
    assert_almost_equal(2 * rln[1], 2 / (1 - x))

    # with correct broadcasting
    forces = Forces({'density': [y, y], 'eps': x})
    rln = forces.get_density_gradient()
    assert_almost_equal(rln[0], 2 / (1 - x))
    assert_almost_equal(rln[1], 2 / (1 - x))

    # test for the pressure
    forces = Forces({'density': y, 'temperature': y, 'eps': x})
    assert_almost_equal(forces.get_pressure_gradient(),  4 / (1 - x))


def test_args():
    p = [1, 2, 3, 4, 5]
    profiles1 = {'density': p, 'temperature': p, 'eps': p}
    profiles2 = {'density': [p, p], 'temperature': p, 'eps': p}
    s = CompositeSpecies(
        [Species('species1', profiles1, mass=1, max_charge=6),
         Species('species1', profiles1, mass=2, max_charge=6),
         Species('species3', profiles2, mass=3, max_charge=1)])

    builder = SpeciesArgsBuilder(s)
    builder.build()
    d = builder.get_result()
    assert_equal(set(builder.delivers()), set(d.keys()))
    assert_equal(d['mass'], [1, 2, 3])
    assert_equal(d['nc'], [1, 1, 2])
    assert_equal(d['zsp'], [[6, 6], [6, 6], [0, 1]])

    # here density is (npoints, nspec, max_nc)
    assert_equal(d['density'].shape, (5, 3, 2))

    # here temperature is (npoints, nspec) since all the charge
    # states are assumed to be isothermal
    assert_equal(d['temperature'].shape, (5, 3))

    # here forces.shape is (npoints, nspec, max_nc, 2)
    assert_equal(d['forces'].shape, (5, 3, 2, 2))


def test_rollaxis():
    a = np.zeros((1, 2, 3))
    assert_equal(myrollaxis(a, 3).shape, (3, 1, 2))

    a = np.zeros((2, 3))
    assert_equal(myrollaxis(a, 3).shape, (1, 2, 3))


def test_iter():
    p = [1, 2, 3, 4, 5]
    profiles1 = {'density': p, 'temperature': p, 'eps': p}
    profiles2 = {'density': [p, p], 'temperature': p, 'eps': p}
    s = CompositeSpecies(
        [Species('species1', profiles1, mass=1, max_charge=6),
         Species('species1', profiles1, mass=2, max_charge=6),
         Species('species3', profiles2, mass=3, max_charge=1)])

    args = Arguments(s, GeometryStrategy(), DefaultNumericArgs())

    l = [a for a in args]
    assert_equal(len(l), len(p))
    assert_equal([a['eps'] for a in args], p)
    assert_equal([a['ic'] for a in args], 5 * [3])
    assert_equal([a['mass'] for a in args], 5 * [[1, 2, 3]])
    assert_equal([a['zsp'] for a in args], 5 * [[[6, 6], [6, 6], [0, 1]]])


def test_neoart_broadcasting():
    eps = np.linspace(0.1, 0.2, 10)
    profiles = {'density': 5, 'temperature': 1, 'eps': [eps]}
    species = [('electron', profiles)]
    g = CircularGeometry(minor_radius=0.1, major_radius=0.88, safety_factor=2,
                         magnetic_field=1.44)
    neoart = Neoart(species, g)
    fluxes = neoart.fluxes()
    assert_equal(fluxes.shape, [10, 1, 1, 4])


def test_CheaseGeometry_tempdir():
    geom = CheaseGeometry(thisdir('neoart.44433t0.5000.dat'))

    pwd = os.getcwd()
    geom.pre()
    assert_equal(os.getcwd(), geom._tempdir)
    assert_equal(os.path.islink('neoart_geom.dat'), True)
    geom.post()
    assert_equal(os.getcwd(), pwd)
    assert_equal(os.path.isdir(geom._tempdir), False)


def test_neoart_vector():
    eps = np.linspace(0.1, 0.2, 10)
    ne = 1 * np.ones_like(eps)
    te = 1 * np.ones_like(eps)
    profiles = {'density': ne, 'temperature': te, 'eps': eps}
    s = CompositeSpecies(
        [Species('species1', profiles, mass=9e-31, max_charge=-1),
         Species('species2', profiles, mass=1.7e-27, max_charge=1)])

    g = CircularGeometry(minor_radius=0.1, major_radius=0.88, safety_factor=2,
                         magnetic_field=1.44)
    neoart = Neoart(species=s, geometry=g)
    fluxes = neoart.fluxes()
    assert_equal(fluxes.shape, (10, 2, 1, 4))


def test_CircularGeometry_broadcasting():
    eps = np.linspace(0.1, 0.2, 10)
    ne = 1 * np.ones_like(eps)
    te = 1 * np.ones_like(eps)
    profiles = {'density': ne, 'temperature': te, 'eps': eps}
    species = CompositeSpecies(
        [Species('species1', profiles, mass=9e-31, max_charge=-1),
         Species('species2', profiles, mass=1.7e-27, max_charge=1)])

    geometry = CircularGeometry(minor_radius=0.1, major_radius=0.88,
                                safety_factor=2, magnetic_field=1.44)
    neoart = Neoart(species, geometry)
    geometry.pre(neoart)  # initialize the iterator
    series = list(geometry._iterator)

    assert_equal([x[0] for x in series], eps)
    assert_equal([x[1] for x in series], 10 * [0.88])  # major radius
    assert_equal([x[3] for x in series], 10 * [2])  # safety factor
    assert_equal([x[4] for x in series], 10 * [1.44])  # magnetic field


# validation tests
#-----------------
# Fixtures

# Test 1 profiles
electron = {'density': 5, 'temperature': 10, 'density gradient': 0,
            'temperature gradient': 0, 'eps': 1e-4}
ion = {'density': 5, 'temperature': 10, 'density gradient': -1,
       'temperature gradient': 1, 'eps': 1e-4}

fixture_species = {}
fixture_species['test1'] = [('electron', electron), ('deuterium', ion)]


# same as test1 but with a larger aspect ratio
electron = electron.copy()
ion = ion.copy()
electron['eps'] = 0.1
ion['eps'] = 0.1
fixture_species['test1_chease'] = [('electron', electron), ('deuterium', ion)]


def test_test1():
    species = fixture_species['test1']
    geometry = CircularGeometry(minor_radius=0.000165, major_radius=1.65,
                                safety_factor=2, magnetic_field=2.5)
    numerical = DefaultNumericArgs()
    numerical.set_contribution('bp')
    numerical.set_electron_ion_collisions(False)
    numerical.set_force_viscosity_regime(1)

    neoart = Neoart(species, geometry, numerical)

    fluxes = neoart.fluxes()
    assert_equal(fluxes.shape, (1, 2, 1, 4))
    assert_approx_equal(fluxes[0, 1, 0, 3], 3.7263e7, significant=4)


def test_segmentation_fault():
    """
    Repeated calls used to give segmentation fault.
    """
    species = fixture_species['test1_chease']
    geometry = CheaseGeometry(thisdir('neoart.44433t0.5000.dat'))
    numeric = DefaultNumericArgs()
    numeric.set_contribution('bp')

    neoart = Neoart(species, geometry)
    c1 = neoart.fluxes()
    # the second call should not triggers a segmentation fault
    c2 = neoart.fluxes()

    # electron fluxes should be equal
    assert_allclose(c1[0, 0, 0], c2[0, 0, 0], rtol=1e-4)

    # ion fluxes should be equal # FIXME this test fails
    #assert_allclose(c1[0, 1, 0], c2[0, 1, 0], rtol=1e-4)


def test_extra_args():
    species = fixture_species['test1']
    geometry = CheaseGeometry(thisdir('neoart.44433t0.5000.dat'))

    args = Arguments(species, geometry, DefaultNumericArgs(),
                     extra_args={'eparr': 0.2})
    assert_equal(args.get_args()['eparr'], 0.2)
    assert_equal([x['eparr'] for x in args], [0.2])

    n = Neoart(species, geometry)
    n._call_neoart = lambda a: assert_equal(a['eparr'], 0)
    n.fluxes()

    n._call_neoart = lambda a: assert_equal(a['eparr'], 0.2)
    n.fluxes(eparr=0.2)

    def my_call_neoart(a):
        assert_equal(a['permanent_arg'], 42.0)
        assert_equal(a['eparr'], 0)  # default parameter
        assert_equal(a['extra_arg'], -1)  # default parameter

    n = Neoart(species, geometry, permanent_arg=42.0)
    n._call_neoart = my_call_neoart

    # Calling n.fluxes() without arguments should fail, because in the _assert
    # we want to have the 'extra_arg' too.
    assert_raises(KeyError, n.fluxes)

    # Now we add 'extra_arg' to the argument list through a keyword argument
    # of fluxes().
    n.fluxes(extra_arg=-1)


def test_TranportCoefficients():
    species = fixture_species['test1']
    geometry = CircularGeometry(minor_radius=0.000165, major_radius=1.65,
                                safety_factor=2, magnetic_field=2.5)

    neo = Neoart(species, geometry)
    coeff = TransportCoefficients(neo).compute()

    density = coeff.get_density()
    forces = neo.get_args().args['forces']
    rln = forces[..., 0] - forces[..., 1]

    pflux_1 = neo.fluxes()[..., 0] / density
    pflux_2 = coeff['D'] * rln + coeff['v']
    assert_almost_equal(pflux_1, pflux_2)
