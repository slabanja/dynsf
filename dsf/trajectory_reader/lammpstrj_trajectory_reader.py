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
import re
from itertools import count
from numpy import array, arange, zeros


class lammpstrj_trajectory_reader(abstract_trajectory_reader):
    """Read LAMMPS trajectory file

    This is a naive (and comparatively slow) implementation,
    written entirely in python.
    """

    @classmethod
    def reader_available(cls):
        return True

    def __init__(self, filename, x_factor=0.1, t_factor=1.0):
        if filename.endswith('.gz'):
            from gzip import GzipFile
            self._fh = GzipFile(filename, 'r')
        elif filename.endswith('.bz2'):
            from bz2 import BZ2File
            self._fh = BZ2File(filename, 'r')
        else:
            self._fh = open(filename, 'r')

        self._open = True
        self._item_re = re.compile(r'^ITEM: (TIMESTEP|NUMBER OF ATOMS|BOX BOUNDS|ATOMS) ?(.*)$')
        self.x_factor = x_factor
        self.t_factor = t_factor
        self.v_factor = x_factor / t_factor
        self._first_called = False
        self._index = count(1)

    # ITEM: TIMESTEP
    # 81000
    # ITEM: NUMBER OF ATOMS
    # 1536
    # ITEM: BOX BOUNDS pp pp pp
    # 1.54223 26.5378
    # 1.54223 26.5378
    # 1.54223 26.5378
    # ITEM: ATOMS id type x y z vx vy vz
    # 247 1 3.69544 2.56202 3.27701 0.00433856 -0.00099307 -0.00486166
    # 249 2 3.73324 3.05962 4.14359 0.00346029 0.00332502 -0.00731005
    # 463 1 3.5465 4.12841 5.34888 0.000523332 0.00145597 -0.00418675

    def _read_frame_header(self):
        while True:
            L = self._fh.readline()
            m = self._item_re.match(L)
            if not m:
                if L == '':
                    self._fh.close()
                    self._open = False
                    raise StopIteration
                if L.strip() == '':
                    continue
                raise IOError("TRJ_reader: Failed to read/parse TRJ frame header")
            if m.group(1) == "TIMESTEP":
                step = int(self._fh.readline())
            elif m.group(1) == "NUMBER OF ATOMS":
                natoms = int(self._fh.readline())
            elif m.group(1) == "BOX BOUNDS":
                bbounds = [map(float, self._fh.readline().split())
                           for _ in range(3)]
                x = array(bbounds)
                box = np.diag(x[:, 1] - x[:, 0])
                if x.shape == (3, 3):
                    box[1, 0] = x[0, 2]
                    box[2, 0] = x[1, 2]
                    box[2, 1] = x[2, 2]
                elif x.shape != (3, 2):
                    raise IOError('TRJ_reader: Malformed box bounds in TRJ frame header')
            elif m.group(1) == "ATOMS":
                cols = tuple(m.group(2).split())
                # At this point, there should be only atomic data left
                return (step, natoms, box, cols)

    def _get_first(self):
        # Read first frame, update state of self, create indexes etc
        step, N, box, cols = self._read_frame_header()
        self._natoms = N
        self._step = step
        self._cols = cols
        self._box = box

        def _all_in_cols(keys):
            for k in keys:
                if not k in cols:
                    return False
            return True

        self._x_map = None
        if _all_in_cols(('id', 'xu', 'yu', 'zu')):
            self._x_I = array(map(cols.index, ('xu', 'yu', 'zu')))
        elif _all_in_cols(('id', 'x', 'y', 'z')):
            self._x_I = array(map(cols.index, ('x', 'y', 'z')))
        elif _all_in_cols(('id', 'xs', 'ys', 'zs')):
            self._x_I = array(map(cols.index, ('xs', 'ys', 'zs')))
            _x_factor = self._box.diagonal().reshape((3, 1))
            # xs.shape == (3,n)
            self._x_map = lambda xs : xs * _x_factor
        else:
            raise RuntimeError('TRJ file must contain at least atom-id, x, y, '
                               'and z coordinates to be useful.')
        self._id_I = cols.index('id')

        if _all_in_cols(('vx', 'vy', 'vz')):
            self._v_I = array(map(cols.index, ('vx', 'vy', 'vz')))
        else:
            self._v_I = None

        if 'type' in cols:
            self._type_I = cols.index('type')
        else:
            self._type_I = None

        data = array([map(float, self._fh.readline().split())
                         for _ in range(N)])
        I = np.asarray(data[:, self._id_I], dtype=np.int)
        # Unless dump is done for group "all" ...
        I[np.argsort(I)] = arange(len(I))
        self._x = zeros((3, N), order='F')
        if self._x_map is None:
            self._x[:, I] = data[:, self._x_I].transpose()
        else:
            self._x[:, I] = self._x_map(data[:, self._x_I].transpose())
        if self._v_I is not None:
            self._v = zeros((3, N), order='F')
            self._v[:, I] = data[:, self._v_I].transpose()

    def _get_next(self):
        # get next frame, update state of self
        step, N, box, cols = self._read_frame_header()
        assert(self._natoms == N)
        assert(self._cols == cols)
        self._step = step
        self._box = box

        data = array([map(float, self._fh.readline().split())
                         for _ in range(N)])
        I = np.asarray(data[:, self._id_I], dtype=np.int) - 1
        if self._x_map is None:
            self._x[:, I] = data[:, self._x_I].transpose()
        else:
            self._x[:, I] = self._x_map(data[:, self._x_I].transpose())
        if self._v_I is not None:
            self._v[:, I] = data[:, self._v_I].transpose()

    def __iter__(self):
        return self

    def close(self):
        if not self._fh.closed:
            self._fh.close()

    def next(self):
        if not self._open:
            raise StopIteration

        if self._first_called:
            self._get_next()
        else:
            self._get_first()

        res = dict(
            index=self._index.next(),
            N=int(self._natoms),
            box=self.x_factor * self._box.copy('F'),
            time=self.t_factor * self._step,
            x=self.x_factor * self._x,
            )

        if self._v_I is not None:
            res['v'] = self.v_factor * self._v
        else:
            res['v'] = None

        return res

