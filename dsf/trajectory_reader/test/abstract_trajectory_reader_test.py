# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.

import numpy
import os

class TrajectoryReaderTestMixin(object):

    LAMMPSTRJ_FIRST_FRAME_FIRST_X = numpy.array([0.191468, 0.302071, 0.0528818])

    def filename_lammpstrj(self):
        data_path = self._data_dir_path()
        return os.path.join(data_path, "positions_and_velocities.lammpstrj")

    def filename_lammpstrj_no_velocities(self):
        data_path = self._data_dir_path()
        return os.path.join(data_path, "positions.lammpstrj")

    def assert_arrays_equal_within_float32eps(self, a, b):
        abs_diff = numpy.absolute(a - b)
        eps = self._float32abs()
        self.assertTrue((abs_diff < eps).all(),
                        "%s and %s not within eps from each other" % (a, b))

    def _data_dir_path(self):
        this_dir = os.path.dirname(__file__)
        return os.path.join(this_dir, "data")

    def _float32abs(self):
        return numpy.finfo(numpy.float32).eps
