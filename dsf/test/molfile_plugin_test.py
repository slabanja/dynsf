import unittest
import dsf.molfile_plugin as mplug


_molfile_not_found = mplug.molfile_plugin_dir() is None

class MolfilePluginTest(unittest.TestCase):
    @unittest.skipIf(_molfile_not_found, "Molfile plugins not found")
    def test_can_load_plugins(self):
        for _, _, _, plugin_name in mplug.TRAJECTORY_PLUGIN_MAPPING:
            plugin_path = mplug.molfile_plugin_path(plugin_name)
            mp = mplug._MolfilePlugin(plugin_path)
            print("%-20s %-15s (%5s %5s %5s %5s)" % \
                      (plugin_name, mp.plugin.filename_extension,
                       bool(mp.plugin.read_timestep),
                       bool(mp.plugin.read_timestep_metadata),
                       bool(mp.plugin.read_next_timestep),
                       bool(mp.plugin.read_structure)))
            mp.close()
