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

import unittest
import dsf.trajectory_reader.molfile_plugin as mplugin


_molfile_not_found = mplugin.molfile_plugin_dir() is None

class MolfilePluginTest(unittest.TestCase):
    @unittest.skipIf(_molfile_not_found, "Molfile plugins not found")
    def test_can_load_plugins(self):
        for _, _, _, plugin_name in mplugin.TRAJECTORY_PLUGIN_MAPPING:
            plugin_path = mplugin.molfile_plugin_path(plugin_name)
            mp = mplugin._MolfilePlugin(plugin_path)
            print("%-20s %-15s (%5s %5s %5s %5s)" % \
                      (plugin_name, mp.plugin.filename_extension,
                       bool(mp.plugin.read_timestep),
                       bool(mp.plugin.read_timestep_metadata),
                       bool(mp.plugin.read_next_timestep),
                       bool(mp.plugin.read_structure)))
            mp.close()
