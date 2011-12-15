import numpy as np


def defaultconfig():
    config = Configuration()
    config.add_species('electron', mass=9.1096e-31, charge=-1)
    config.add_species('ion', mass=1.6727e-27, charge=1)

    config.set_temperature(10)
    config.set_density(5, ['ion', 'electron'])

    return config


class Configuration(object):
    def __init__(self):
        self.species = {}
        self._species_names = []
        self.temperatures = {}
        self.densities = {}

        self.pressure_gradients = {}
        self.temperature_gradients = {}

        self.loop_voltage = 0.
        self.accuracy = 1e-5
        self.calculate_geometry = True
        self.contribution = 3

        self.numerical_parameters = {
            'nreg' : 1,
            'nleg' : 3,
            'nenergy' : 1,
            'ncof' : 0,
            'neofrc' : 0,
            'sigma' : np.zeros(4),
            'ishot' : 0,
            'isel' : 2,
        }

    def add_species(self, name, mass, charge):
        self._species_names.append(name)
        self.species[name] = Species(name, mass, charge)

    def set_temperature(self, temperature, species='all'):
        """
        Set the temperature of a species in [keV].  If species is 'all' then
        all the species are set accordingly.
        """
        for s in self._parse_argument(species):
            self.temperatures[s] = temperature

    def set_density(self, density, species='all'):
        """
        Set the density of a species in [10^19 m^-3].  If species is 'all'
        then all the species are set accordingly.
        """
        for s in self._parse_argument(species):
            self.densities[s] = density

    def set_pressure_gradient(self, pprim, species='all'):
        for s in self._parse_argument(species):
            self.pressure_gradients[s] = pprim

    def set_temperature_gradient(self, tprim, species='all'):
        for s in self._parse_argument(species):
            self.temperature_gradients[s] = tprim

    def set_contribution(self, contribution='all'):
        c = contribution.lower()
        if c in ['classical', 'c']:
            self.contribution = 0
        elif c in ['banana-plateau', 'bp']:
            self.contribution = 1
        elif c in ['pfirsch-schlueter', 'ps']:
            self.contribution = 2
        elif c == 'all':
            self.contribution = 3
        else:
            raise ValueError('invalid transport regime.')

    def _parse_argument(self, species):
        if species == 'all': return self._species_names
        if isinstance(species, str): return [species]
        return species

    def todict(self):
        num_of_species = len(self.species)
        num_of_ionisation_stages = np.ones(num_of_species)

        zsp = np.zeros((num_of_species, 1))

        temperatures = np.zeros(num_of_species)
        den = np.zeros((num_of_species, 1))
        masses = np.zeros(num_of_species)
        ds = np.zeros((num_of_species,1, 2))

        for i, s in enumerate(self._species_names):
            masses[i] = self.species[s].mass
            zsp[i, 0] = self.species[s].charge
            temperatures[i] = self.temperatures[s]
            den[i, 0] = self.densities[s]

            ds[i, 0, 0] = self.pressure_gradients.get(s, 0.)
            ds[i, 0, 1] = self.temperature_gradients.get(s, 0.)

        d = dict(nc=num_of_ionisation_stages, zsp=zsp, masses=masses,
            den=den, temperatures=temperatures, ds=ds)

        d.update(self.numerical_parameters)

        d['eparr'] = self.loop_voltage
        d['eps'] = self.accuracy
        d['neogeo'] = int(self.calculate_geometry)
        d['ic'] = self.contribution
        return d


class Species(object):
    def __init__(self, name, mass, charge):
        self.name = name
        self.mass = mass
        self.charge = charge

