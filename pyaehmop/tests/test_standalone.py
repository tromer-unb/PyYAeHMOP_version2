"""Tests use Python only; never call YAeHMOP C."""
import unittest
from pathlib import Path
import numpy as np
from pyaehmop import read_input, calculate
from pyaehmop.basis import Orbital
from pyaehmop.overlap import overlap_block
from pyaehmop.engine import occupancies

DATA = Path(__file__).resolve().parents[1] / "inputs"


class TestStandalone(unittest.TestCase):
    def test_analytic_slater_1s_overlap(self):
        zeta=1.3
        orb=(Orbital(0,"s",1,0,-13.6,(zeta,),(1.,)),)
        r=0.74*1.889644746
        expected=np.exp(-zeta*r)*(1+zeta*r+(zeta*r)**2/3)
        measured=overlap_block(orb,orb,[0,0,0.74])[0,0]
        self.assertAlmostEqual(measured,expected,places=9)

    def test_hermiticity_e_symmetry(self):
        s=read_input(DATA / "XeF4")
        calc=calculate(s)
        S=calc["overlap_matrix"]
        H=calc["hamiltonian_matrix"]
        self.assertTrue(np.allclose(S,S.T,atol=1e-10))
        self.assertTrue(np.allclose(H,H.T,atol=1e-10))
        self.assertTrue(np.allclose(calc["coefficients"].T@S@calc["coefficients"],
                                    np.eye(calc["orbitals"]),atol=1e-9))
        self.assertTrue(np.all(np.linalg.eigvalsh(S)>0))

    def test_xef4_public_c_reference(self):
        calc=calculate(read_input(DATA / "XeF4"))
        self.assertEqual(calc["orbitals"],20)
        self.assertAlmostEqual(calc["electronic_energy_eV"],-812.312883,places=3)
        self.assertAlmostEqual(calc["homo_eV"],-11.344305,places=3)
        self.assertAlmostEqual(calc["lumo_eV"],0.835790,places=3)
        self.assertAlmostEqual(calc["gap_homo_lumo_eV"],12.18009,places=3)
        self.assertAlmostEqual(np.sum(calc["mulliken_charges"]),0.,places=7)
        self.assertAlmostEqual(calc["mulliken_charges"][0],3.40524,places=3)

    def test_feh_double_zeta_d(self):
        calc=calculate(read_input(DATA / "FeH"))
        self.assertEqual(calc["orbitals"],10)
        self.assertAlmostEqual(calc["electronic_energy_eV"],-118.442593,places=3)
        self.assertAlmostEqual(calc["gap_homo_lumo_eV"],0.6232,places=3)
        self.assertAlmostEqual(sum(calc["mulliken_charges"]),0.,places=6)

    def test_diamond_serial_parallel_bands(self):
        s=read_input(DATA / "diamond")
        serial=calculate(s,workers=1)
        parallel=calculate(s,workers=2)
        self.assertEqual(serial["orbitals"],8)
        self.assertEqual(len(serial["kpoints"]),81)
        self.assertTrue(np.allclose(serial["band_energies_eV"],
                                    parallel["band_energies_eV"],atol=1e-9))
        self.assertAlmostEqual(serial["band_path_gap_eV"],4.9605249,delta=0.001)

    def test_coop_molecular(self):
        calc=calculate(read_input(DATA / "XeF4"))
        curves=calc["bond_analyses"]
        self.assertEqual(len(curves),4)
        for item in curves:
            self.assertAlmostEqual(item["integrated_coop"],0.014274,places=5)

    def test_periodic_properties_and_dos(self):
        calc=calculate(read_input(DATA / "diamond"),workers=2)
        av=calc["sampled_properties"]
        self.assertAlmostEqual(av["average_electronic_energy_eV"],-134.924903,places=3)
        self.assertAlmostEqual(av["fermi_eV"],-9.370505,places=3)
        self.assertAlmostEqual(av["charge_sum"],0.,places=7)
        self.assertAlmostEqual(av["dos"]["integrated_states"],16,places=4)

    def test_degeneracy(self):
        occ=occupancies([-5.,-5.,0.,0.],3)
        np.testing.assert_allclose(occ,[1.5,1.5,0.,0.])


if __name__ == "__main__":
    unittest.main()
