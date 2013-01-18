from collections import Iterable, Sequence
import numpy as np


smart_species = False
try:
    import periodictable
    import scipy.constants as constants
    smart_species = True
except ImportError:
    raise Warning('automatic species identification is not possible')


class ProfileCollection(object):
    """
    container for different profiles: density, temperature, rotation etc.

    always converts to ndarray
    """

    def __init__(self, data=None, broadcast=True):
        self._profiles = {}

        if isinstance(data, dict):
            self._init_dict(data)
        elif isinstance(data, self.__class__):
            self._init_CompositeProfile(data)

        if broadcast:
            self.broadcast()

    def _init_dict(self, data):
        for key, value in data.iteritems():
            self.add(key, value)

    def _init_CompositeProfile(self, data):
        self._profiles = data._profiles

    def add(self, name, p):
        self._profiles[name] = np.asarray(p)

    def remove(self):
        pass

    def get_profile(self, name):
        return self._profiles.get(name, np.array([]))

    def broadcast(self):
        if self._profiles == {}:
            return
        broadcasted = np.broadcast_arrays(*self._profiles.values())
        self._profiles = dict(zip(self._profiles.keys(), broadcasted))

    def has_profile(self, name):
        return name in self._profiles.keys()


class MultistageProfileCollection(ProfileCollection):

    def get_profile(self, name):
        """
        Always returns a profile with at least 2 dimensions, where the first
        dimension represents the number of ionisation stages.
        """
        p = self._profiles.get(name, np.array([]))
        p = np.atleast_1d(p)
        if (p.ndim == 1):
            p = p[np.newaxis]
        return p


# Composite pattern: SpeciesComponent, Species, CompositeSpeices
class SpeciesComponent(object):
    """
    Abstract class to define the SpeciesComponent interface.
    """

    def get_mass(self):
        return self._mass

    def get_max_charge(self):
        return self._max_charge

    def get_name(self):
        return self._name

    def get_profile(self):
        """
        return the certain profile of a certain ionisation stage.  If None and
        the profile.ndim == 1 it is assumed that this belongs to the fully
        ionized stage.

        0 neutral, None fully ionized or integer
        """
        pass

    mass = property(get_mass)
    max_charge = property(get_max_charge)


class Species(SpeciesComponent):

    def __init__(self, name, profiles=None, mass=None, max_charge=None):
        if isinstance(name, Species):
            mass = name.get_mass()
            max_charge = name.get_max_charge()
            profiles = name._profiles
            name = name.name
        elif isinstance(name, Sequence) and isinstance(name[1], dict):
            name, profiles = name

        self._name = name
        self._profiles = MultistageProfileCollection(profiles)
        if (mass is None) or (max_charge is None):
            mass = self.guess_mass()
            max_charge = self.guess_max_charge()
        self._mass = mass
        self._max_charge = max_charge

    def get_profile(self, name):
        return self._profiles.get_profile(name)

    def get_zsp(self, name):
        p = self.get_profile(name)
        # if there's only one ionisation stage it is considered as the fully
        # ionized one
        if (p.shape[0] == 1):
            return np.array([self.max_charge])
        elif (p.shape[0] == abs(self.max_charge) + 1):
            return np.arange(abs(self.max_charge) + 1,
                             step=np.sign(self.max_charge))
        else:
            raise ValueError('zsp')

    def guess_mass(self):
        if self.name in ['electron']:
            return constants.electron_mass
        else:
            element = self.lookup_element()
            if element is None:
                return None
            else:
                return element.mass * constants.atomic_mass

    def guess_max_charge(self):
        if self.name in ['electron']:
            return -1
        else:
            element = self.lookup_element()
            if element is None:
                return None
            else:
                return element.number

    def has_profile(self, name):
        return self._profiles.has_profile(name)

    def lookup_element(self):
        if not smart_species:
            return None
        return getattr(periodictable, self.name, None)

    profiles = property(lambda self: self._profiles)
    name = property(lambda self: self._name)


class CompositeSpecies(SpeciesComponent):

    def __init__(self, data=None):
        self._species = []

        if isinstance(data, CompositeSpecies):
            self._species = data._species
        elif isinstance(data, Iterable):
            self._init_iterable(data)

    def _init_iterable(self, data):
        for s in data:
            self.add(s)

    def add(self, s):
        if isinstance(s, Sequence):
            s = Species(name=s[0], profiles=s[1])
        self._species.append(s)

    def remove(self, name):
        pass

    def get_max_charge(self):
        return [s.get_max_charge() for s in self._species]

    def get_mass(self):
        return [s.get_mass() for s in self._species]

    def get_zsp(self, name):
        return [s.get_zsp(name) for s in self._species]

    def get_profile(self, name):
        return [s.get_profile(name) for s in self._species]

    def __iter__(self):
        return iter(self._species)

    def get_names(self):
        return [s.get_name() for s in self._species]

    names = property(get_names)
