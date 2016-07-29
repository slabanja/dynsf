
# Copyright (C) 2011 Mattias Slabanja <slabanja@chalmers.se>
#
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

"""Molfile plugin interface

Not much more than the necessary ctypes stuff.
"""

import os
import re
from sys import maxsize
from platform import uname
from itertools import islice

from ctypes import (
    cdll, CDLL, RTLD_GLOBAL, POINTER, CFUNCTYPE,
    Structure, cast, pointer, c_int, c_uint, c_float, c_double, c_char, c_char_p, c_void_p)
from ctypes.util import find_library

_system = uname()[0]
_64bit_python = maxsize > 2 ** 32
if _system != "Windows":
    # a bit ugly. Is there a more kosher way of making sure that the
    # libstdc++ (or similar) symbols are available to the molfile-plugins?
    _cxx = CDLL(find_library('stdc++'), mode=RTLD_GLOBAL)


c_int_p = POINTER(c_int)
c_int_pp = POINTER(c_int_p)
c_char_pp = POINTER(c_char_p)
c_char_ppp = POINTER(c_char_pp)
c_float_p = POINTER(c_float)
c_float_pp = POINTER(c_float_p)

#
# As of molfile_plugin abiversion 16, the following trajectory
# formats are supported (pasted from the AMD web page):
#
# Molecular Dynamics Trajectory File Plugins
#
#     AMBER 'binpos' trajectory reader (.binpos)
#     AMBER "CRD" trajectory reader (.crd, .crdbox)
#     CHARMM, NAMD, X-PLOR "DCD" reader/writer (.dcd)
#     CPMD (CPMD trajectory) reader (.cpmd)
#     DLPOLY HISTORY file reader (.dlpolyhist)
#     Gromacs TRR/XTC reader (.trr, .xtc)
#     LAMMPS trajectory reader (.lammpstrj)
#     MMTK NetCDF trajectory reader (.nc)
#     VASP trajectories of ionic steps (.xml, .OUTCAR, .XCATCAR)
#     VTF trajectory files (.vtf)
#     XCrySDen, Quantum Espresso XSF/AXSF trajectory files (.axsf, .xsf)
#     XYZ trajectory files (.xyz)
#

TRAJECTORY_PLUGIN_MAPPING = (
    # (SOFTWARE-NAME, FILE-TYPE, FILE-SUFFIX, PLUGIN-NAME)
    ('AMBER', '"binpos"',
     'binpos', 'binposplugin'),

    ('AMBER', '"CRD"',
     'crd', 'crdplugin'),

    ('CHARMM', '"DCD" - CHARMM, NAMD, XPLOR',
     'dcd', 'dcdplugin'),

    ('CPMD', 'CPMD',
     'cpmd', 'cpmdplugin'),

    ('DLPOLY', 'DLPOLY History',
     'dlpolyhist', 'dlpolyplugin'),

    ('GROMACS', 'Gromacs XTC',
     'xtc', 'gromacsplugin'),

    ('GROMACS', 'Gromacs TRR',
     'trr', 'gromacsplugin'),

    ('LAMMPS', 'LAMMPS Trajectory',
     'lammpstrj', 'lammpsplugin'),

    ('VASP', 'VASP ionic steps',
     'xml', 'vaspxmlplugin'),

    ('VASP', 'VASP ionic steps',
     'OUTCAR', 'vaspoutcarplugin'),

    ('VTF', 'VTF trajectory',
     'vtf', 'vtfplugin'),

    ('XCrySDen', 'XSF trajectory',
      'xsf', 'xsfplugin'),

    ('XCrySDen', 'AXSF trajectory',
      'axsf', 'xsfplugin'),

    ('?', 'XYZ trajectory',
     'xyz', 'xyzplugin'))


MIN_ABI_VERSION = 15

MOLFILE_PLUGIN_TYPE = "mol file reader"
VMDPLUGIN_SUCCESS = 0
VMDPLUGIN_ERROR = -1
class vmdplugin_t(Structure):
    _fields_ = [('abiversion', c_int),
                ('type', c_char_p),
                ('name', c_char_p),
                ('prettyname', c_char_p),
                ('author', c_char_p),
                ('majorv', c_int),
                ('minorv', c_int),
                ('is_reentrant', c_int)]

class molfile_metadata_t(Structure):
    _fields_ = [('database', c_char * 81),
                ('accession', c_char * 81),
                ('date', c_char * 81),
                ('title', c_char * 81),
                ('remarklen', c_int),
                ('remarks', c_char_p)]

class molfile_atom_t(Structure):
    _fields_ = [('name', c_char * 16),
                ('type', c_char * 16),
                ('resname', c_char * 16),
                ('resid', c_int),
                ('segid', c_char * 16),
                ('chain', c_char * 16),
                ('altloc', c_char * 16),
                ('insertion', c_char * 16),
                ('occupancy', c_float),
                ('bfactor', c_float),
                ('mass', c_float),
                ('charge', c_float),
                ('radius', c_float),
                ('atomicnumber', c_int)]

class molfile_timestep_metadata_t(Structure):
    _fields_ = [('count', c_uint),
                ('avg_bytes_per_timestamp', c_uint),
                ('has_velocities', c_int)]

class molfile_qm_metadata_t(Structure):
    # Left as an exercise to the reader
    pass

class molfile_qm_timestep_t(Structure):
    # Left as an exercise to the reader
    pass

class molfile_timestep_t(Structure):
    _fields_ = [('coords', POINTER(c_float)),
                ('velocities', POINTER(c_float)),
                ('A', c_float),
                ('B', c_float),
                ('C', c_float),
                ('alpha', c_float),
                ('beta', c_float),
                ('gamma', c_float),
                ('physical_time', c_double)]

class molfile_volumetric_t(Structure):
    _fields_ = [('dataname', c_char * 256),
                ('origin', c_float * 3),
                ('xaxis', c_float * 3),
                ('yaxis', c_float * 3),
                ('zaxis', c_float * 3),
                ('xsize', c_int),
                ('ysize', c_int),
                ('zsize', c_int),
                ('has_color', c_int)]


dummy_fun_t = CFUNCTYPE(c_int)
class molfile_plugin_t(Structure):
    # ABI from molfile abiversion 16
    # This is the important(TM) structure.
    # Only partial read support is considered for now, other
    # functions have been replaced by a dummy_fun_t placeholder.
    _fields_ = [('abiversion', c_int),
                ('type', c_char_p),
                ('name', c_char_p),
                ('prettyname', c_char_p),
                ('author', c_char_p),
                ('majorv', c_int),
                ('minorv', c_int),
                ('is_reentrant', c_int),
                ('filename_extension', c_char_p),
#
# void *(* open_file_read)(const char *filepath, const char *filetype, int *natoms);
                ('open_file_read', CFUNCTYPE(c_void_p, c_char_p, c_char_p, c_int_p)),
#
# int (*read_structure)(void *, int *optflags, molfile_atom_t *atoms);
                ('read_structure', CFUNCTYPE(c_int, c_void_p, c_int_p,
                                       POINTER(molfile_atom_t))),
#
# int (*read_bonds)(void *, int *nbonds, int **from, int **to, float **bondorder,
#                   int **bondtype, int *nbondtypes, char ***bondtypename);
                ('read_bonds', dummy_fun_t),
#
# int (* read_next_timestep)(void *, int natoms, molfile_timestep_t *);
                ('read_next_timestep', CFUNCTYPE(c_int, c_void_p, c_int,
                                           POINTER(molfile_timestep_t))),
#
# void (* close_file_read)(void *);
                ('close_file_read', CFUNCTYPE(None, c_void_p)),
#
# void *(* open_file_write)(const char *filepath, const char *filetype,
#      int natoms);
                ('open_file_write', dummy_fun_t),
#
#  int (* write_structure)(void *, int optflags, const molfile_atom_t *atoms);
                ('write_structure', dummy_fun_t),
#
#  int (* write_timestep)(void *, const molfile_timestep_t *);
                ('write_timestep', dummy_fun_t),
#
#  void (* close_file_write)(void *);
                ('close_file_write', dummy_fun_t),
#
#  int (* read_volumetric_metadata)(void *, int *nsets,
#        molfile_volumetric_t **metadata);
                ('read_volumetric_metadata', CFUNCTYPE(c_int, c_void_p, c_int_p,
                                                 POINTER(POINTER(molfile_volumetric_t)))),
#
#  int (* read_volumetric_data)(void *, int set, float *datablock,
#        float *colorblock);
                ('read_volumetric_data', CFUNCTYPE(c_int, c_void_p, c_int, c_float_p,
                                             c_float_p)),
#
#  int (* read_rawgraphics)(void *, int *nelem, const molfile_graphics_t **data);
                ('read_rawgraphics', dummy_fun_t),
#
#  int (* read_molecule_metadata)(void *, molfile_metadata_t **metadata);
                ('read_molecule_metadata', CFUNCTYPE(c_int, c_void_p,
                                               POINTER(POINTER(molfile_metadata_t)))),
#
#  int (* write_bonds)(void *, int nbonds, int *from, int *to, float *bondorder,
#                     int *bondtype, int nbondtypes, char **bondtypename);
                ('write_bonds', dummy_fun_t),
#
#  int (* write_volumetric_data)(void *, molfile_volumetric_t *metadata,
#                                float *datablock, float *colorblock);
                ('write_volumetric_data', dummy_fun_t),
#
#  int (* read_angles)(void *handle, int *numangles, int **angles, int **angletypes,
#                      int *numangletypes, char ***angletypenames, int *numdihedrals,
#                      int **dihedrals, int **dihedraltypes, int *numdihedraltypes,
#                      char ***dihedraltypenames, int *numimpropers, int **impropers,
#                      int **impropertypes, int *numimpropertypes, char ***impropertypenames,
#                      int *numcterms, int **cterms, int *ctermcols, int *ctermrows);
                ('read_angles', dummy_fun_t),
#
#  int (* write_angles)(void *handle, int numangles, const int *angles, const int *angletypes,
#                       int numangletypes, const char **angletypenames, int numdihedrals,
#                       const int *dihedrals, const int *dihedraltypes, int numdihedraltypes,
#                       const char **dihedraltypenames, int numimpropers,
#                       const int *impropers, const int *impropertypes, int numimpropertypes,
#                       const char **impropertypenames, int numcterms, const int *cterms,
#                       int ctermcols, int ctermrows);
                ('write_angles', dummy_fun_t),
#
#  int (* read_qm_metadata)(void *, molfile_qm_metadata_t *metadata);
                ('read_qm_metadata', dummy_fun_t),
#
#  int (* read_qm_rundata)(void *, molfile_qm_t *qmdata);
                ('read_qm_rundata', dummy_fun_t),
#
#  int (* read_timestep)(void *, int natoms, molfile_timestep_t *,
#                        molfile_qm_metadata_t *, molfile_qm_timestep_t *);
                ('read_timestep', CFUNCTYPE(c_int, c_void_p, c_int, POINTER(molfile_timestep_t),
                                      POINTER(molfile_qm_metadata_t),
                                      POINTER(molfile_qm_timestep_t))),
#
#  int (* read_timestep_metadata)(void *, molfile_timestep_metadata_t *);
                ('read_timestep_metadata', CFUNCTYPE(c_int, c_void_p,
                                               POINTER(molfile_timestep_metadata_t))),
#
#  int (* read_qm_timestep_metadata)(void *, molfile_qm_timestep_metadata_t *);
                ('read_qm_timestep_metadata', dummy_fun_t),
#
#  int (* cons_fputs)(const int, const char*);
                ('cons_fputs', CFUNCTYPE(c_int, c_int, c_char_p))]




# typedef int (*vmdplugin_register_cb)(void *, vmdplugin_t *);
vmdplugin_register_cb_t = CFUNCTYPE(c_int, c_void_p, POINTER(vmdplugin_t))


class _MolfilePlugin(object):
    """A thin molfile_plugin wrapper class

    This class holds the loaded plugin-library and sets up
    a molfile_plugin_t structure.

    The 'close' method, calls the plugin fini-function.
    """
    def __init__(self, plugin_path):
        library = self._load_plugin_library(plugin_path)

        if library.vmdplugin_init() != VMDPLUGIN_SUCCESS:
            raise RuntimeError('Failed to init %s' % plugin_path)

        plugin_p = pointer(molfile_plugin_t(0))
        callback = self._create_vmdplugin_register_cb(plugin_p)

        if library.vmdplugin_register(None, callback) != VMDPLUGIN_SUCCESS:
            raise RuntimeError('Failed to register %s' % plugin_path)

        self._library = library
        self.plugin = plugin_p.contents

    def close(self):
        self.plugin = molfile_plugin_t(0)
        self._library.vmdplugin_fini()

    def _load_plugin_library(self, plugin_path):
        library = cdll.LoadLibrary(plugin_path)

        # extern int vmdplugin_init(void);
        library.vmdplugin_init.restype = c_int

        # extern int vmdplugin_fini(void);
        library.vmdplugin_fini.restype = c_int

        # extern int vmdplugin_register(void *, vmdplugin_register_cb);
        library.vmdplugin_register.restype = c_int
        library.vmdplugin_register.argtypes = (c_void_p, vmdplugin_register_cb_t)

        return library

    def _create_vmdplugin_register_cb(self, plugin_p):
        def py_register_callback(v, p):
            contents = p.contents
            if contents.type == MOLFILE_PLUGIN_TYPE:
                if contents.abiversion >= MIN_ABI_VERSION:
                    plugin_p.contents = cast(p, POINTER(molfile_plugin_t)).contents
            return VMDPLUGIN_SUCCESS

        vmdplugin_register_cb = vmdplugin_register_cb_t(py_register_callback)

        return vmdplugin_register_cb


def molfile_plugin_path(plugin_name):
    plugin_dir = molfile_plugin_dir()
    plugin_suffix = _molfile_machine_specific_plugin_suffix()
    path = os.path.join(plugin_dir, plugin_name + "." + plugin_suffix)
    if os.path.exists(path):
        return path

def molfile_plugin_dir():
    vmddir = _vmddir()
    if vmddir is not None:
        return os.path.join(vmddir, _molfile_relative_plugin_dir())

def _molfile_relative_plugin_dir():
    "Molfile path relative VMD root"
    return os.path.join('plugins',
                        _molfile_machine_specific_plugin_dir_name(),
                        'molfile')

def _molfile_machine_specific_plugin_suffix():
    return 'so'

def _molfile_machine_specific_plugin_dir_name():
    if _system == 'Linux':
        if _64bit_python:
            return 'LINUXAMD64'
        else:
            return 'LINUX'

    elif not _64bit_python:
        if _system == 'Darwin':
            return 'MACOSXX86'

        elif _system == 'Windows':
            return 'WIN32'

    bits = "64" if _64bit_python else "32"
    raise RuntimeError("No known %s-bit plugin-path on %s system" % (bits, _system))


def _vmddir():
    vmddir = os.environ.get('VMDDIR')
    if vmddir is not None:
        return vmddir

    for path in os.environ['PATH'].split(os.pathsep):
        tentative_vmd_script_path = os.path.join(path, 'vmd')
        if os.path.exists(tentative_vmd_script_path):
            vmddir = _get_vmddir_from_vmd_script(tentative_vmd_script_path)
            if vmddir:
                return vmddir

    for tentive_vmddir in ["C:\Program Files (x86)\University of Illinois\VMD"]:
        if os.path.exists(tentive_vmddir):
            return tentive_vmddir


def _get_vmddir_from_vmd_script(vmd_script_path):
    with open(vmd_script_path) as fh:
        for L in islice(fh, 10):
            match = re.match(r'^(?:set )?defaultvmddir=(?:"(/.*)"|(/[^# ]*)).*$', L)
            if match:
                a, b = match.groups()
                return a or b
