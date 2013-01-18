"""
This module computes transport coefficients.
"""
import numpy as np

from arguments import SpeciesArgsBuilder


class TransportCoefficients(object):
    def __init__(self, neoart):
        self.neoart = neoart

        self.coeffs = {}
        for key in ['D', 'v']:
            self.coeffs[key] = np.zeros_like(self.get_density())

    def compute(self):
        nc = self.neoart.get_args().args['nc']
        D = self.coeffs['D']
        v = self.coeffs['v']

        forces = self.neoart.get_args().args['forces']
        density = self.get_density()

        for isp, s in enumerate(nc):
            for ic in xrange(s):
                pflux = self.flux_with_density_gradient_only(isp, ic)
                D[:, isp, ic] = (pflux / density)[:, isp, ic]

                # contribution to all other fluxes appears as pinch
                for jsp, r in enumerate(nc):
                    for jc in xrange(r):
                        if (jc != ic) or (isp != jsp):
                            rln = forces[:, isp, ic, 0] - forces[:, isp, ic, 1]
                            v[:, jsp, jc] += (pflux[:, jsp, jc]
                                              * rln / density[:, jsp, jc])

            v += self.flux_with_temperature_gradient_only(isp, ic) / density

        return self

    def flux_with_density_gradient_only(self, isp, ic):
        """
        Compute particle flux where only the isp, ic density gradient is
        not zero.
        """
        args = self.neoart.get_args()
        args.builders['species'] = OnlyDensityGradient(self.neoart._species,
                                                       isp, ic)
        args.build()

        return self.particle_flux(args)

    def flux_with_temperature_gradient_only(self, isp, ic):
        args = self.neoart.get_args()
        args.builders['species'] = \
            OnlyTemperatureGradient(self.neoart._species, isp, ic)
        args.build()

        return self.particle_flux(args)

    def particle_flux(self, args, **kwargs):
        fluxes = []
        self.neoart._geometry.pre(self.neoart)
        for a in args:
            self.neoart._geometry.at_each_rho(a)
            fluxes.append(self.neoart._call_neoart(a))
        self.neoart._geometry.post(self.neoart)

        return np.array(fluxes)[..., 0]  # particle flux

    def get_density(self):
        p = np.array(self.neoart._species.get_profile('density'))
        return np.rollaxis(p, -1) * 1e19

    def __getitem__(self, key):
        return self.coeffs[key]


class OnlyOneSpeciesGradient(SpeciesArgsBuilder):
    def __init__(self, species, isp, ic):
        super(OnlyOneSpeciesGradient, self).__init__(species)
        self.isp = isp
        self.ic = ic


class OnlyDensityGradient(OnlyOneSpeciesGradient):
    def _build_forces(self):
        d = super(OnlyDensityGradient, self)._build_forces()
        d['forces'].fill(0)
        d['forces'][:, self.isp, self.ic, 0] = 1
        return d


class OnlyTemperatureGradient(OnlyOneSpeciesGradient):
    def _build_forces(self):
        d = super(OnlyTemperatureGradient, self)._build_forces()
        forces = np.zeros_like(d['forces'])
        forces[:, self.isp, self.ic, 0] = d['forces'][:, self.isp, self.ic, 1]
        forces[:, self.isp, self.ic, 1] = d['forces'][:, self.isp, self.ic, 1]
        d['forces'] = forces
        return d
