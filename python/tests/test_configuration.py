import numpy as np
import unittest

from neoart.configuration import Configuration


class TestConfiguration(unittest.TestCase):
    def testAddSpecies(self):
        config = Configuration()
        config.add_species('electron', mass=9.1096e-31, charge=-1)
        config.add_species('ion', mass=1.6727e-27, charge=1)

        self.assertEqual(len(config.species), 2)
        self.assertEqual(config._species_names[0], 'electron')
        self.assertEqual(config.species['electron'].mass, 9.1096e-31)
        self.assertEqual(config.species['electron'].charge, -1)
        self.assertEqual(config._species_names[1], 'ion')


    def testSetTemperature(self):
        config = Configuration()
        config.add_species('electron', mass=1, charge=-1)
        config.add_species('ion', mass=1, charge=1)

        config.set_temperature(1, 'electron')
        self.assertEqual(config.temperatures['electron'], 1)

        config.set_temperature(2, ['ion', 'electron'])
        self.assertEqual(config.temperatures['electron'], 2)
        self.assertEqual(config.temperatures['ion'], 2)
        
        config.set_temperature(3)
        self.assertEqual(config.temperatures['electron'], 3)
        self.assertEqual(config.temperatures['ion'], 3)

    def testArgumentParsing(self):
        config = Configuration()
        config.add_species('electron', mass=1, charge=-1)
        config.add_species('ion', mass=1, charge=1)

        d = config._parse_argument('electron')
        self.assertEqual(d, ['electron'])

        d = config._parse_argument('all')
        self.assertEqual(set(d), set(['electron', 'ion']))

        d = config._parse_argument(['ion', 'electron'])
        self.assertEqual(set(d), set(['electron', 'ion']))

    def testSetPressureGradient(self):
        config = Configuration()
        config.add_species('electron', mass=1, charge=-1)
        config.add_species('ion', mass=1, charge=1)

        config.set_pressure_gradient(123, 'ion')
        self.assertEqual(config.pressure_gradients['ion'], 123)
        self.assertEqual(config.pressure_gradients.get('electron', None),
                None)

    def testTodict(self):
        config = Configuration()
        config.add_species('electron', mass=1, charge=-1)
        config.add_species('ion', mass=2, charge=1)

        config.set_temperature(10)
        config.set_density(5)

        d = config.todict()
        np.testing.assert_array_equal(d['masses'], [1., 2.])
        np.testing.assert_array_equal(d['temperatures'], [10., 10.])
        np.testing.assert_array_equal(d['den'], [[5.], [5.]])
        np.testing.assert_array_equal(d['nc'], [1, 1])


if __name__ == '__main__':
    unittest.main()
