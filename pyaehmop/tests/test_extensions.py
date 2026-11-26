"""Extension regressions: Z-matrix, f, double-zeta, FMO, COHP and exports.

Never runs YAeHMOP C; reference numbers were obtained once
from the original BSD-licensed distribution.
"""
import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np

from pyaehmop import read_input, calculate
from pyaehmop.basis import angular_functions
from pyaehmop.geometry import zmat_to_cartesian
from pyaehmop.export import export_csv

DATA = Path(__file__).resolve().parents[1] / "inputs"


class TestExtensions(unittest.TestCase):
    def test_water_zmatrix_against_original(self):
        s=read_input(DATA/"water_zmat")
        self.assertEqual(len(s.atoms),3)
        np.testing.assert_allclose(s.atoms[1].position,[0,0,.9572],atol=1e-12)
        np.testing.assert_allclose(s.atoms[2].position,[.926627206486,0,-.239987208409],
                                   atol=1e-10)
        report=calculate(s)
        self.assertAlmostEqual(report["electronic_energy_eV"],-162.541767,places=3)
        self.assertAlmostEqual(report["lumo_eV"],-.188029,places=3)

    def test_zmatrix_documented_dihedral(self):
        geometria=zmat_to_cartesian([
            "1 Bi","2 C 1 2.1","3 I 1 2.7 2 90",
            "4 I 1 2.7 2 90 3 90","5 I 1 2.7 2 90 3 180"])
        pos=[x[1] for x in geometria]
        np.testing.assert_allclose(pos[3],[0,-2.7,0],atol=1e-12)
        np.testing.assert_allclose(pos[4],[-2.7,0,0],atol=1e-12)

    def test_ceh_4f(self):
        res=calculate(read_input(DATA/"CeH"))
        self.assertEqual(res["orbitals"],17)
        self.assertAlmostEqual(res["electronic_energy_eV"],-62.191213,places=3)
        self.assertAlmostEqual(res["homo_eV"],-11.28,places=5)
        self.assertAlmostEqual(sum(res["mulliken_charges"]),0,places=7)

    def test_smh_4f_double_zeta(self):
        s=read_input(DATA/"SmH")
        self.assertEqual(len([x for x in s.atoms[0].shells if x.l==3]),1)
        f=next(x for x in s.atoms[0].shells if x.l==3)
        self.assertEqual(len(f.zetas),2)
        res=calculate(s)
        self.assertAlmostEqual(res["electronic_energy_eV"],-107.40948,places=3)
        self.assertAlmostEqual(res["lumo_eV"],-11.2475,places=3)

    def test_orthonormal_f_harmonics(self):
        mu,w=np.polynomial.legendre.leggauss(24)
        phi=np.arange(64)*2*np.pi/64
        sin=np.sqrt(np.maximum(0,1-mu**2))
        xyz=np.stack(np.broadcast_arrays(sin[:,None]*np.cos(phi)[None,:],
                                         sin[:,None]*np.sin(phi)[None,:],
                                         mu[:,None]+0*phi[None,:]),axis=-1)
        Y=np.stack(angular_functions(xyz,3),axis=-1)
        gram=np.einsum("ipm,ipn,i->mn",Y,Y,w)*2*np.pi/len(phi)
        np.testing.assert_allclose(gram,np.eye(7),atol=2e-12)

    def test_xef4_fmo(self):
        res=calculate(read_input(DATA/"XeF4"))
        fmo=res["fmo"]
        self.assertEqual(len(fmo["fragments"]),2)
        self.assertEqual([x["orbital_count"] for x in fmo["fragments"]],[4,16])
        self.assertAlmostEqual(fmo["fragments"][0]["fragment_electronic_energy_eV"],
                               -66.68,places=5)
        self.assertAlmostEqual(fmo["fragments"][1]["fragment_electronic_energy_eV"],
                               -754.357,places=3)
        S=fmo["fmo_overlap_matrix"]
        np.testing.assert_allclose(np.diag(S),np.ones(20),atol=1e-10)

    def test_coop_cohp_periodic(self):
        res=calculate(read_input(DATA/"diamond_coop"),workers=2)
        group=res["bond_groups"][0]
        self.assertAlmostEqual(group["average_integrated_coop"],0.702514,places=5)
        self.assertAlmostEqual(group["average_integrated_cohp"],-17.016028,places=4)

    def test_bonds_intercell(self):
        res=calculate(read_input(DATA/"diamond_intercell"),workers=2)
        bonds=res["bond_analyses"]
        self.assertEqual(len(bonds),3)
        self.assertAlmostEqual(bonds[0]["integrated_coop"],-.000798,places=5)
        self.assertAlmostEqual(bonds[1]["integrated_coop"],.515354,places=5)
        self.assertAlmostEqual(bonds[2]["integrated_cohp"],.045501,places=5)

    def test_parallelism_molecular(self):
        s=read_input(DATA/"XeF4")
        r1=calculate(s,workers=1)
        r2=calculate(s,workers=2)
        np.testing.assert_allclose(r1["orbital_energies_eV"],r2["orbital_energies_eV"],atol=1e-10)
        np.testing.assert_allclose(r1["mulliken_charges"],r2["mulliken_charges"],
                                   atol=1e-10)

    def test_csv_export(self):
        mol=calculate(read_input(DATA/"water_zmat"))
        solid=calculate(read_input(DATA/"diamond"),workers=1)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            files=export_csv(mol,root/"molecule")
            self.assertEqual(len(files),3)  # levels, Mulliken charges and energies
            with (root/"molecule"/"molecular_levels.csv").open() as fd:
                rows=list(csv.reader(fd))
            self.assertEqual(len(rows),mol["orbitals"]+1)
            paths=export_csv(solid,root/"solid")
            self.assertEqual(len(paths),4)  # bands, DOS, Mulliken e energias
            with (root/"solid"/"bands.csv").open() as fd:
                rows=list(csv.reader(fd))
            self.assertEqual(len(rows),1+81*8)


if __name__=="__main__":
    unittest.main()
