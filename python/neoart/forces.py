import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

from species import MultistageProfileCollection


class Forces(object):
    """
    Thermodynamic forces
    """

    def __init__(self, profiles, gradient_strategy=None, squeeze=True):
        self.profiles = MultistageProfileCollection(profiles)
        self.set_gradient_strategy(gradient_strategy)
        self.squeeze = squeeze

    def get_pressure_gradient(self):
        return self.get_temperature_gradient() + self.get_density_gradient()

    def get_temperature_gradient(self):
        return self.scale_length('temperature')

    def get_density_gradient(self):
        return self.scale_length('density')

    def set_gradient_strategy(self, strategy):
        if strategy is None:
            strategy = DefaultStrategy(self.profiles)
        self.gradient_strategy = strategy

    def scale_length(self, name):
        ret = self.gradient_strategy.scale_length(name)
        if self.squeeze:
            ret = ret.squeeze()
        return ret


# Strategies to calculate the profile gradients
class GradientStrategy(object):

    def __init__(self, profiles):
        """
        Parameters
        ----------
        profile : Species object
        """
        self.profiles = profiles

    def scale_length(self, name):
        pass

    def applicable(self, name):
        """
        Return True if the strategy is applicable.
        """
        return True


class CopyStrategy(GradientStrategy):

    def source_label(self, name):
        return name + ' gradient'

    def scale_length(self, name):
        return self.profiles.get_profile(self.source_label(name))

    def applicable(self, name):
        return self.profiles.has_profile(self.source_label(name))


class InterpolatingStrategy(GradientStrategy):

    def __init__(self, profiles, xlabel='eps'):
        super(InterpolatingStrategy, self).__init__(profiles)
        self.xlabel = xlabel

    def scale_length(self, name):
        """
        compute - 1/y dydx
        """
        y = self.profiles.get_profile(name)
        x = self.profiles.get_profile(self.xlabel)

        if x.shape != y.shape:
            return np.array([])

        res = np.zeros_like(x)
        for i in xrange(x.shape[0]):
            yfit = InterpolatedUnivariateSpline(x[i], y[i])
            res[i] = - yfit(x[i], 1) / yfit(x[i])

        return res

    def applicable(self, name):
        return self.profiles.has_profile(self.xlabel)


class DefaultStrategy(GradientStrategy):

    def __init__(self, profiles):
        super(DefaultStrategy, self).__init__(profiles)
        self.copy = CopyStrategy(profiles)
        self.interpolate = InterpolatingStrategy(profiles)

    def scale_length(self, name):
        if self.copy.applicable(name):
            return self.copy.scale_length(name)
        elif self.interpolate.applicable(name):
            return self.interpolate.scale_length(name)
        else:
            raise Warning("huhhh, gradients?")
