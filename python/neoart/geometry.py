import os
import tempfile
import numpy as np

from itertools import izip

import core


class GeometryStrategy(object):
    """
    Abstract class to define the geometry interface.
    """

    def get_args(self):
        """
        return the neoart input parameters of a particular geometry
        """
        return {}

    def pre(self, neoart=None):
        """
        command to invoke before running neoart
        """
        pass

    def post(self, neoart=None):
        """
        command to invoke after running neoart
        """
        pass

    def at_each_rho(self, args=None):
        pass


class CircularGeometry(GeometryStrategy):
    def __init__(self, minor_radius, major_radius, safety_factor,
                 magnetic_field):
        """
        amin : minor_radius
        rmag : major_radius
        q : safety factor
        bn : magnetic field strength
        """

        self.minor_radius = minor_radius
        self.major_radius = major_radius
        self.safety_factor = safety_factor
        self.magnetic_field = magnetic_field
        self.inverse_aspect_ratio = self.minor_radius / self.major_radius

        self._iterator = None

    def at_each_rho(self, args=None):
        params = self._iterator.next()
        core.circgeom(*params)

    def pre(self, neoart):
        eps = np.array([x['eps'] for x in neoart.get_args()])
        broadcasted = np.broadcast_arrays(eps, self.major_radius,
                                          self.inverse_aspect_ratio,
                                          self.safety_factor,
                                          self.magnetic_field)
        self._iterator = izip(*broadcasted)

    def get_args(self):
        return {'isel': 2, 'rho': 0}


class CheaseGeometry(GeometryStrategy):

    NEOART_GEOM = 'neoart_geom.dat'

    def __init__(self, chease_file):
        if not os.path.isfile(chease_file):
            raise IOError("No such file: '%s'" % chease_file)
        self.chease_file = os.path.abspath(chease_file)

    def pre(self, neoart=None):
        self._change_to_tempdir()
        self._link_chease_file_to_tempdir()

    def post(self, neoart=None):
        self._restore_pwd()

    def get_args(self):
        return {'isel': 4}

    def _change_to_tempdir(self):
        self._cwd = os.getcwd()
        self._tempdir = tempfile.mkdtemp(suffix='.neoart')
        os.chdir(self._tempdir)

    def _restore_pwd(self):
        os.chdir(self._cwd)
        os.remove(os.path.join(self._tempdir, self.NEOART_GEOM))
        os.rmdir(self._tempdir)

    def _link_chease_file_to_tempdir(self):
        os.symlink(self.chease_file,
                   os.path.join(self._tempdir, self.NEOART_GEOM))
