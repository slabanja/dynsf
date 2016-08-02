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

from ctypes import cdll, POINTER, c_float, c_int, c_char_p
from ctypes.util import find_library
from itertools import count
from dsf.trajectory_reader.abstract_trajectory_reader import abstract_trajectory_reader
import numpy as np

#
# L I B G M X
#
# libgmx comes with Gromacs and contains, among other things,
# functionality for reading xtc-files.
#

libgmx_name = find_library('gromacs')
libgmx = libgmx_name and cdll.LoadLibrary(libgmx_name)
np_ndp = np.ctypeslib.ndpointer

if libgmx:
    # single prec gmx-real equals float, right?
    xtcfloat_np = np.float32
    xtcfloat_ct = c_float
    xtcint_ct = c_int

    # t_fileio *open_xtc(const char *filename,const char *mode);
    # /* Open a file for xdr I/O */
    libgmx.open_xtc.restype = POINTER(xtcint_ct)
    libgmx.open_xtc.argtypes = [c_char_p, c_char_p]

    # int read_first_xtc(t_fileio *fio,
    #                           int *natoms,int *step,real *time,
    #                           matrix box,rvec **x,real *prec,gmx_bool *bOK);
    # /* Open xtc file, read xtc file first time, allocate memory for x */
    libgmx.read_first_xtc.restype = xtcint_ct
    libgmx.read_first_xtc.argtypes = [
        POINTER(xtcint_ct), POINTER(xtcint_ct),
        POINTER(xtcint_ct), POINTER(xtcfloat_ct),
        np_ndp(dtype=xtcfloat_np, shape=(3, 3),
               flags='f_contiguous, aligned'),
        POINTER(POINTER(xtcfloat_ct)),
        POINTER(xtcfloat_ct), POINTER(xtcint_ct)]

    # int read_next_xtc(t_fileio *fio,
    #                          int natoms,int *step,real *time,
    #                          matrix box,rvec *x,real *prec,gmx_bool *bOK);
    # /* Read subsequent frames */
    libgmx.read_next_xtc.restype = xtcint_ct
    libgmx.read_next_xtc.argtypes = [
        POINTER(xtcint_ct), xtcint_ct,
        POINTER(xtcint_ct), POINTER(xtcfloat_ct),
        np_ndp(dtype=xtcfloat_np, shape=(3, 3),
               flags='f_contiguous, aligned'),
        np_ndp(dtype=xtcfloat_np, ndim=2,
               flags='f_contiguous, aligned'),
        POINTER(xtcfloat_ct), POINTER(xtcint_ct)]


class xtc_trajectory_reader(abstract_trajectory_reader):

    @classmethod
    def reader_available(cls):
        # TODO: Fix so that this works also for gromacs 5
        # When changing from libgmx to libgromacs-5,
        # it seems the library function signatures now
        # include an XDR-data structore.
        # Investigate, implement, and test...
        return False
#         return libgmx is not None

    def __init__(self, filename):
        if libgmx is None:
            raise RuntimeError("XTC_reader: No libgmx found, can't use XTC_reader!")

        self._fio = libgmx.open_xtc(filename, 'rb')
        if not self._fio:
            raise IOError("XTC_reader: Failed to open file %s" % filename)

        self._index = count(1)
        self._natoms = xtcint_ct()
        self._step = xtcint_ct()
        self._time = xtcfloat_ct()
        self._box = np.require(np.zeros((3, 3)),
                                  xtcfloat_np, ['F_CONTIGUOUS', 'ALIGNED'])
        self._x = None
        self._prec = xtcfloat_ct()
        self._bOK = xtcint_ct()  # gmx_bool equals int
        self._open = True
        self._first_called = False

    def _get_first(self):
        # Read first frame and update state of self accordingly
        _xfirst = POINTER(xtcfloat_ct)()
        res = libgmx.read_first_xtc(self._fio, self._natoms,
                                    self._step, self._time,
                                    self._box, _xfirst,
                                    self._prec, self._bOK)
        self._first_called = True
        if not res:
            raise IOError("XTC_reader: read_first_xtc failed")
        if not self._bOK.value:
            raise IOError("XTC_reader: corrupt frame in xtc-file?")

        N = self._natoms.value
        self._x = np.require(array(_xfirst[0:3 * N]).reshape((3, N), order='F'),
                             xtcfloat_np, ['F_CONTIGUOUS', 'ALIGNED'])
        self._x.flags.writeable = False

    def _get_next(self):
        # get next frame, update state of self
        res = libgmx.read_next_xtc(self._fio, self._natoms.value,
                                   self._step, self._time,
                                   self._box, self._x,
                                   self._prec, self._bOK)
        if not res:
            return False
        if not self._bOK.value:
            raise IOError("XTC_reader: corrupt frame in xtc-file?")
        return True

    def __iter__(self):
        return self

    def close(self):
        if self._open:
            libgmx.close_xtc(self._fio)
            self._open = False

    def next(self):
        if not self._open:
            raise StopIteration

        if self._first_called:
            if not self._get_next():
                self.close()
                raise StopIteration
        else:
            self._get_first()

        return dict(
            index=self._index.next(),
            box=self._box.copy('F'),
            time=self._time.value,
            N=self._natoms.value,
            x=self._x,
            v=None,
            )

