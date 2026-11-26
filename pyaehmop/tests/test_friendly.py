"""Tests for human-friendly XYZ/CIF/ASE/YAeHMOP input, without invoking C."""
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np

from pyaehmop.friendly import (
    read_easy_job,calculate_easy_job,is_easy_input
)
from pyaehmop.input import read_input
from pyaehmop.engine import calculate
from pyaehmop.export import export_csv

ROOT=Path(__file__).resolve().parents[1]
EXAMPLES=ROOT/"examples_friendly"


class TestFriendlyInput(unittest.TestCase):
    def _job(self,name):
        return read_easy_job(EXAMPLES/(name+".job"))

    def test_auto_detection(self):
        self.assertTrue(is_easy_input(EXAMPLES/"water_charges.job"))
        self.assertTrue(is_easy_input(EXAMPLES/"diamond_cif.job"))
        self.assertFalse(is_easy_input(ROOT/"inputs"/"water_zmat"))
        self.assertFalse(is_easy_input(ROOT/"inputs"/"XeF4"))

    def test_xyz_only_charges(self):
        report=calculate_easy_job(self._job("water_charges"))
        self.assertEqual(report["general_input"]["requested_calculations"],["charges"])
        self.assertNotIn("orbital_energies_eV",report)
        self.assertNotIn("dos",report)
        self.assertNotIn("electronic_energy_eV",report)
        self.assertNotIn("band_energies_eV",report)
        self.assertAlmostEqual(float(sum(report["mulliken_charges"])),0,places=10)
        with tempfile.TemporaryDirectory() as folder:
            files=export_csv(report,folder)
            self.assertEqual({Path(x).name for x in files},{"mulliken.csv"})

    def test_xyz_matches_zmatrix_energy(self):
        report=calculate_easy_job(self._job("water_full"))
        c=calculate(read_input(ROOT/"inputs"/"water_zmat"))
        self.assertAlmostEqual(report["electronic_energy_eV"],
                               c["electronic_energy_eV"],places=7)
        self.assertAlmostEqual(report["homo_eV"],-14.8,places=4)
        self.assertIn("dos",report)
        self.assertEqual(len(report["dos"]["energy_eV"]),401)
        self.assertAlmostEqual(report["dos"]["sigma_eV"],0.2)
        with tempfile.TemporaryDirectory() as folder:
            files=export_csv(report,folder)
            self.assertEqual({Path(x).name for x in files},
                             {"molecular_levels.csv","dos.csv",
                              "mulliken.csv","energies.csv"})

    def test_yaehmop_xef4_includes_fmo_and_coop(self):
        report=calculate_easy_job(self._job("xef4_original"))
        self.assertEqual(report["general_input"]["source_format"],"yaehmop")
        self.assertAlmostEqual(report["electronic_energy_eV"],-812.312883,places=3)
        self.assertEqual(len(report["bond_analyses"]),4)
        self.assertEqual(len(report["fmo"]["fragments"]),2)

    def test_yaehmop_diamond_unchanged(self):
        report=calculate_easy_job(self._job("diamond_reference"))
        self.assertEqual(report["type"],"periodic")
        self.assertEqual(len(report["kpoints"]),81)
        props=report["sampled_properties"]
        self.assertAlmostEqual(props["average_electronic_energy_eV"],
                               -134.924903,places=4)
        self.assertAlmostEqual(props["charge_sum"],0,places=10)

    @unittest.skipUnless(__import__("importlib").util.find_spec("ase"),
                         "Optional ASE package is not installed")
    def test_cif_matches_custom_yaehmop(self):
        reference=calculate_easy_job(self._job("diamond_reference"))
        cif=calculate_easy_job(self._job("diamond_cif"))
        self.assertEqual(len(cif["band_energies_eV"]),81)
        np.testing.assert_allclose(cif["band_energies_eV"],
                                   reference["band_energies_eV"],atol=1e-6)
        self.assertAlmostEqual(
            cif["sampled_properties"]["average_electronic_energy_eV"],
            reference["sampled_properties"]["average_electronic_energy_eV"],
            places=6)

    @unittest.skipUnless(__import__("importlib").util.find_spec("ase"),
                         "Optional ASE package is not installed")
    def test_periodic_extxyz_charges_energy(self):
        report=calculate_easy_job(self._job("diamond_ase"))
        self.assertNotIn("band_energies_eV",report)
        self.assertNotIn("dos",report)
        props=report["sampled_properties"]
        self.assertEqual(props["total_kpoints"],8)
        self.assertEqual(len(props["mulliken_charges"]),2)
        self.assertLess(abs(sum(props["mulliken_charges"])),1e-8)

    @unittest.skipUnless(__import__("importlib").util.find_spec("ase"),
                         "Optional ASE package is not installed")
    def test_automatic_ase_bandpath(self):
        source=EXAMPLES/"diamond_primitive.cif"
        with tempfile.TemporaryDirectory() as folder:
            jobfile=Path(folder)/"auto.job"
            jobfile.write_text(f"cif\n{source}\ncalculations = bands\n"
                               "band_path = auto\ndivisions = 35\n"
                               "images = 2 2 2\n")
            job=read_easy_job(jobfile)
            self.assertEqual(job.source_type,"cif")
            self.assertEqual(len(job.bands_kpoints),35)
            self.assertGreater(len(job.bands_labels),1)
            self.assertFalse(job.structure.kpoints)

    def test_clear_core_error_message(self):
        source=EXAMPLES/"water.xyz"
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/"invalid.job"
            p.write_text(f"xyz\n{source}\ncalculations=bands\n")
            with self.assertRaisesRegex(ValueError,"periodic structure"):
                read_easy_job(p)
            p.write_text(f"xyz\n{source}\nfoo=bar\n")
            with self.assertRaisesRegex(ValueError,"unknown option"):
                read_easy_job(p)
            p.write_text(f"xyz\n{source}\ncalculations=charges\n"
                         "charge=0\nelectrons=8\n")
            with self.assertRaisesRegex(ValueError,"OR"):
                read_easy_job(p)

    def test_fragments_in_xyz(self):
        source=EXAMPLES/"water.xyz"
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/"fmo.job"
            p.write_text(f"xyz\n{source}\ncalculations=fmo\n"
                         "fragments=1:6;2-3:2\n")
            report=calculate_easy_job(read_easy_job(p))
            self.assertEqual(len(report["fmo"]["fragments"]),2)
            self.assertNotIn("electronic_energy_eV",report)

    def test_xyz_molecular_bonds(self):
        source=EXAMPLES/"water.xyz"
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/"bonds.job"
            p.write_text(f"xyz\n{source}\ncalculations=bonds\n"
                         "bonds=1-2; 1-3\n")
            report=calculate_easy_job(read_easy_job(p))
            self.assertEqual(len(report["bond_analyses"]),2)
            c=report["bond_analyses"]
            self.assertAlmostEqual(c[0]["integrated_coop"],
                                   c[1]["integrated_coop"],places=7)


    def test_config_explicit_format_and_filename(self):
        source=EXAMPLES/"water.xyz"
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"chaves.job"
            path.write_text(f"format = xyz\nfile = {source}\n"
                            "calculations = only_charges\n"
                            "plots = yes\noutputs = json\n"
                            "workers = 2\n")
            job=read_easy_job(path)
            self.assertEqual(job.source_type,"xyz")
            self.assertIn("charges",job.tasks)
            self.assertIn("png",job.outputs)

    @unittest.skipUnless(__import__("importlib").util.find_spec("ase"),
                         "Optional ASE package is not installed")
    def test_vasp_poscar_via_ase(self):
        job=self._job("diamond_poscar")
        self.assertEqual(job.source_type,"ase")
        self.assertEqual(len(job.structure.atoms),2)
        self.assertEqual(job.structure.bounds,(5,5,5))
        rep=calculate_easy_job(job)
        self.assertNotIn("band_energies_eV",rep)
        self.assertEqual(len(rep["sampled_properties"]["mulliken_charges"]),2)

    @unittest.skipUnless(__import__("importlib").util.find_spec("ase"),
                         "Optional ASE package is not installed")
    def test_periodic_only_dos_without_bands(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"dos.job"
            source=EXAMPLES/"diamond_primitive.cif"
            path.write_text(f"cif\n{source}\ncalculations = dos\n"
                            "kmesh = 2 2 2\ncutoff = 8\n"
                            "images = 4 4 4\noutputs = json,csv\n")
            rep=calculate_easy_job(read_easy_job(path))
            self.assertNotIn("band_energies_eV",rep)
            self.assertNotIn("mulliken_charges",rep.get("sampled_properties",{}))
            self.assertIn("dos",rep["sampled_properties"])
            with tempfile.TemporaryDirectory() as destination:
                files=export_csv(rep,destination)
                self.assertEqual({Path(f).name for f in files},{"dos.csv"})

    @unittest.skipUnless(__import__("importlib").util.find_spec("matplotlib"),
                         "Optional Matplotlib package is not installed")
    def test_output_png(self):
        from pyaehmop.plots import export_plots
        report=calculate_easy_job(self._job("water_full"))
        with tempfile.TemporaryDirectory() as folder:
            paths=export_plots(report,folder)
            self.assertEqual({Path(f).name for f in paths},{"dos.png","charges.png"})
            for file in paths:
                self.assertEqual(Path(file).read_bytes()[:8],b"\x89PNG\r\n\x1a\n")

    def test_yaehmop_format_preserves_data(self):
        source=ROOT/"inputs"/"XeF4"
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"default.job"
            path.write_text("yaehmop\n"+str(source)+"\n")
            job=read_easy_job(path)
            self.assertIn("fmo",job.tasks)
            self.assertIn("bonds",job.tasks)
            self.assertEqual(len(job.structure.atoms),5)

if __name__=="__main__":
    unittest.main()
