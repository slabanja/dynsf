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

import numpy as np
from dsf.trajectory_reader.abstract_trajectory_reader import abstract_trajectory_reader
from dsf.trajectory_reader.molfile_plugin import (
    _MolfilePlugin, molfile_plugin_dir, molfile_plugin_path, TRAJECTORY_PLUGIN_MAPPING,
    molfile_atom_t, molfile_timestep_metadata_t, molfile_timestep_t)
from ctypes import c_float, c_int, byref, POINTER
from itertools import count
from numpy import zeros, pi, array, cos


# Molfile plugin is using single precission floats
molfile_float_np = np.float32
molfile_float_ct = c_float

class molfile_trajectory_reader(abstract_trajectory_reader):
    """Read a trajectory using the molfile_plugin package

    molfile_plugin is a part of VMD, and consists of
    plugins for a fairly large number of different trajectory
    formats (see molfile_plugin.TRAJECTORY_PLUGIN_MAPPING).

    filename - string, filename of trajectory file.
    index_file - string, filename of ini-style index file.
    plugin - string, name of plugin to use. If None, guess
             pluginname by looking at filename suffix.
    """

    @classmethod
    def reader_available(cls):
        return molfile_plugin_dir() is not None

    def _guess_plugin(self, filename):
        filename_suffix = filename.rsplit('.', 1)[-1]
        for _, _, suffix, plugin_name in TRAJECTORY_PLUGIN_MAPPING:
            if suffix == filename_suffix:
                return plugin_name

    def __init__(self, filename, plugin_name=None, x_factor=0.1, t_factor=1.0):

        if plugin_name is None:
            plugin_name = self._guess_plugin(filename)

        if plugin_name is None:
            raise RuntimeError('molfile_reader: no suitable plugin known for file %s' % filename)

        self.x_factor = x_factor
        self.t_factor = t_factor
        self.v_factor = x_factor / t_factor

        self._N = c_int()
        suffix = filename.rsplit('.', 1)[-1]

        self._mfp = _MolfilePlugin(molfile_plugin_path(plugin_name))
        p = self._mfp.plugin

        self._fh = p.open_file_read(filename, suffix,
                                    byref(self._N))
        if not self._fh:
            raise RuntimeError('molfile_reader: failed to open file %s with plugin %s.' % (
                    filename, plugin_name))
        N = self._N.value

        if p.read_structure:
            # for e.g. lammpsplugin, read_structure needs to be called first so
            # that molfile_reader knows (internally) which coordinates and
            # velocites to map to which atoms.
            self._atoms_arr = (molfile_atom_t * N)()
            self._optflags = c_int()
            rc = p.read_structure(self._fh, byref(self._optflags), self._atoms_arr)
            if rc:
                raise IOError('molfile_reader: read structure failed for '
                              'file %s (plugin %s, rc %i)' % (filename, plugin_name, rc))
        else:
            self._atoms_arr = None

        self._v = None
        self._x = np.require(zeros((3, N)), molfile_float_np,
                             ['F_CONTIGUOUS', 'ALIGNED'])
        self._x.flags.writeable = False

        if p.read_timestep_metadata:
            # It seems only lammpsplugin offers this (but other formats could
            # include velocity information. how to test for that!?)
            tsm = self._timestep_metadata = molfile_timestep_metadata_t()
            rc = p.read_timestep_metadata(self._fh, byref(tsm))
            if rc:
                raise IOError('molfile_reader: read timestep metadata failed for '
                              'file %s (plugin %s, rc %i)' % (filename, plugin_name, rc))

            if tsm.has_velocities:
                self._v = np.require(zeros((3, N)), molfile_float_np,
                                     ['F_CONTIGUOUS', 'ALIGNED'])
                self._v.flags.writeable = False
        else:
            self._timestep_metadata = None

        # Now, set up the timestep structure
        self._ts = molfile_timestep_t()
        self._ts.coords = self._x.ctypes.data_as(POINTER(molfile_float_ct))
        if self._v is not None:
            self._ts.velocities = self._v.ctypes.data_as(POINTER(molfile_float_ct))
        else:
            # Set velocities to a NULL pointer
            self._ts.velocities = POINTER(molfile_float_ct)()

        # Set frame counter
        self._index = count(1)

    def __iter__(self):
        return self

    def next(self):
        if not self._mfp.plugin.read_next_timestep:
            raise StopIteration

        ts = self._ts
        rc = self._mfp.plugin.read_next_timestep(self._fh, self._N, byref(ts))
        if rc:
            self._mfp.close()
            raise StopIteration

        res = dict(
                   index=self._index.next(),
                   box=self._to_box(ts.A, ts.B, ts.C,
                                ts.alpha, ts.beta, ts.gamma) * self.x_factor,
                   N=self._N.value,
                   time=ts.physical_time * self.t_factor,
                   x=self._x * self.x_factor
                   )
        if self._v is not None:
            res['v'] = self._v * self.v_factor
        else:
            res['v'] = None

        return res

    def close(self):
        self._mfp.close()


    def _to_box(self, A, B, C, alpha, beta, gamma):
        # create box vectors out of molfile-info
        deg2rad = pi / 180.0
        return array(((A, 0.0, 0.0),
                      (B * cos(deg2rad * gamma), B, 0.0),
                      (C * cos(deg2rad * beta), C * cos(deg2rad * alpha), C)))
