import unittest
import numpy as np

from neoart.neoart import neoart

class TestNeoart(unittest.TestCase):
    """Low level wrapper."""
    def setUp(self):
        nc = [1, 1]
        zsp = [[-1], [1]]
        mass = [9.1096e-31, 1.6727e-27]
        density = [[5], [5]]
        temperature = [10, 10]

        forces = np.zeros((2,1,2))
        forces[1,0,1] = 1 # ion temperature gradient

        g = dict(rho=0.000165, Rmag=1.65, q=2, bn=2.5)
        self.coeff = neoart(nc, zsp, mass, temperature, density, forces,
                geometry=g, contribution='bp', electron_ion_collisions=False,
                force_viscosity=1)

    def test1(self):
        np.testing.assert_approx_equal(self.coeff[1, 0, 3], 3.7263e7,
                significant=4)


from neoart import colxi
class TestCollisionFrequency(unittest.TestCase):
    def test(self):
        nc = [1, 1]
        zsp = [[-1], [1]]
        mass = [9.1096e-31, 1.6727e-27]
        density = [[5], [5]]
        temperature = [10, 10]

        col, xi = colxi(nc, zsp, mass, density, temperature)

        col_desired = np.array([[  1.15759798e-07,   1.17738693e-07],
                                [  5.04520228e-06,   5.96420518e-06]])
        np.testing.assert_array_almost_equal(col, col_desired)


import math
from neoart.neoart import neoart2, circgeom, colxi2
from neoart.configuration import Configuration

test1 = Configuration()
test1.add_species('electron', mass=9.1096e-31, charge=-1.0)
test1.add_species('ion', mass=1.6727e-27, charge=1.0)
test1.set_density(5.) # in 10^19 m^-3
test1.set_temperature(10.) # in keV
test1.set_contribution('BP')
 
rho = 0.05
eps = 1e-4
q = 2.
Rmag = 1.65
bn = 2.5

class TestNeoart2(unittest.TestCase):
    def setUp(self):
        circgeom(rho=rho, Rmag=Rmag, eps=eps, q=q, bn=bn)

        test1.set_pressure_gradient(0); test1.set_temperature_gradient(0)
        test1.set_pressure_gradient(1, 'electron')
        self.cff1 = neoart2(rho, test1)

        test1.set_pressure_gradient(0); test1.set_temperature_gradient(0)
        test1.set_pressure_gradient(1, 'ion')
        self.cff2 = neoart2(rho, test1)

        test1.set_pressure_gradient(0); test1.set_temperature_gradient(0)
        test1.set_temperature_gradient(1, 'electron')
        self.cff3 = neoart2(rho, test1)

        test1.set_pressure_gradient(0); test1.set_temperature_gradient(0)
        test1.set_temperature_gradient(1, 'ion')
        self.cff4 = neoart2(rho, test1)

        self.tau, xi = colxi2(test1)
        self.T = test1.temperatures['ion']

    def test_ion_heat_flux(self):
        """
        Test #1 (test1.f) from the original NEOART source.
        """
        norm = 1.6e-22 * bn**2 * math.sqrt(eps**3)
        norm /= (2 * q**2 * self.T * self.tau[1,1])
        sq2 = math.sqrt(2)

        ion_temperature_gradient = sq2 * self.cff4.energy_flux['ion'] * norm
        electron_temperature_gradient = sq2 * self.cff3.energy_flux['ion'] * norm
        ion_pressure_gradient = sq2 * self.cff2.energy_flux['ion'] * norm
        electron_pressure_gradient = sq2 * self.cff1.energy_flux['ion'] * norm

        self.assertClose(tau[1,1], 5.964e-6)
        self.assertClose(norm, 2.0958e-24)
        self.assertClose(ion_temperature_gradient, -0.68)
        self.assertEqual(electron_temperature_gradient, 0)
        self.assertEqual(ion_pressure_gradient, -2.06930e-15)
        self.assertEqual(electron_pressure_gradient, 0)

    def test_ion_heat_flux(self):
        pass

    def assertClose(self, a, b, accuracy=1e-5):
        self.assertTrue(1 - a/b < accuracy)


if __name__ == '__main__':
    unittest.main()
