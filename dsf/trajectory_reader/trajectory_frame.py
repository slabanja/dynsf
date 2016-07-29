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

class trajectory_frame(object):
    """Trivial data struct holding MD-data for one time frame
    
     'index' : trajectory index,
     'box'   : simulation box as 3 row vectors (nm),
     'N'     : number of atoms,
     'x'     : particle positions as 3xN array (nm),
     'v'     : (*) particle velocities as 3xN array (nm/ps),
     'time'  : (*) simulation time (ps),
     
     (*) may not be available, depends on reader and trajectory file format.
    """
    def __init__(self, **kwargs):
        for attribute, value in kwargs.iteritems():
            setattr(self, attribute, value)

