"""Regression checks for origin-free device CIFs and legacy input compatibility."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import numpy as np
from pyaehmop.device_job import read_device_cif,run_job
from pyaehmop.parameter_files import load_parameter_selection
try:
    import ase
except ImportError:
    ase=None

@unittest.skipIf(ase is None,"ASE is not installed")
class DeviceCIFOriginTest(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.params=load_parameter_selection("default",self.root)

    def cif(self,name,z,length=10):
        p=self.root/name
        p.write_text(f"""data_region
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
_cell_length_a 10
_cell_length_b 10
_cell_length_c {length}
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
H1 H 0.5 0.5 {z/length} 1
""")
        return p

    def test_no_metadata_and_legacy_offset(self):
        p=self.cif("left.cif",2)
        a=read_device_cif(p,self.params)
        b=read_device_cif(p,self.params,{"cartesian_origin_A":[-5,-5,-3],"box_lengths_A":[10,10,10]})
        np.testing.assert_allclose(a.atoms[0].position,[5,5,2],atol=1e-12)
        np.testing.assert_allclose(b.atoms[0].position,[0,0,-1],atol=1e-12)

    def test_mismatched_cells_rejected_before_transport(self):
        self.cif("left.cif",2)
        self.cif("device.cif",4)
        self.cif("right.cif",6,12)
        p=self.root/"test.device.job"
        p.write_text("""structure_format = cif
electrode_left = left.cif
scattering_region = device.cif
electrode_right = right.cif
left_repeat = 0 0 2
right_repeat = 0 0 2
""")
        with patch("pyaehmop.device_job.build_device") as solver:
            with self.assertRaisesRegex(ValueError,"same cell"):
                run_job(p)
            solver.assert_not_called()
