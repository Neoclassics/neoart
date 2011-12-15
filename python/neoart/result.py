class Result(object):
    def __init__(self, out_array, config):
        self.config = config
        self.particle_flux = {}
        self.energy_flux = {}
        self.parallel_flow = {}
        self.poloidal_flow = {}

        for i, s in enumerate(config._species_names):
            self.particle_flux[s] = out_array[i, 0]
            self.energy_flux[s] = out_array[i, 1]
            self.parallel_flow[s] = out_array[i, 2]
            self.poloidal_flow[s] = out_array[i, 3]

    def get_species(self):
        return self.config._species_names

    species = property(get_species)

