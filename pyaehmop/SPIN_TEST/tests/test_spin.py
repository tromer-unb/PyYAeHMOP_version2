"""Independent spin-SCF/SOC tests, without any C executable."""
from dataclasses import replace
from functools import lru_cache
from pathlib import Path
import tempfile
import unittest
import numpy as np

from pyaehmop.basis import Orbital
from pyaehmop.SPIN_TEST.config import read_spin_job
from pyaehmop.SPIN_TEST.operators import angular_momentum,soc_matrix
from pyaehmop.SPIN_TEST.solver import (
    fermi_dirac,spinor_problem,scf,evaluate_bands,spin_dos)
from pyaehmop.SPIN_TEST.plotting import comparison_plot, plot_soc_zoom

ROOT=Path(__file__).resolve().parents[1]


@lru_cache(maxsize=4)
def case(label):
    job=read_spin_job(ROOT/"examples"/f"Fe_wire_{label}.spin")
    res,ref=scf(job.structure,job.kpoints,job.weights,job.settings,
                quadrature=job.quadrature,cutoff=job.cutoff)
    bands,polar=evaluate_bands(job.structure,job.band_path,ref)
    return job,res,ref,bands,polar


class TestOperators(unittest.TestCase):
    def test_angular_momentum_spdf(self):
        for l in range(4):
            lx,ly,lz=angular_momentum(l)
            n=2*l+1
            np.testing.assert_allclose(lx,lx.conj().T,atol=1e-10)
            np.testing.assert_allclose(ly,ly.conj().T,atol=1e-10)
            np.testing.assert_allclose(lz,lz.conj().T,atol=1e-10)
            np.testing.assert_allclose(lx@ly-ly@lx,1j*lz,atol=1e-10)
            np.testing.assert_allclose(lx@lx+ly@ly+lz@lz,
                                       l*(l+1)*np.eye(n),atol=1e-10)

    def test_analytic_p_soc_multiplets(self):
        p=[Orbital(0,label,2,1,-5,(1.,),(1.,))
           for label in ("px","py","pz")]
        xi=.12
        h=soc_matrix(p,["C"],{"C":{"p":xi}})
        np.testing.assert_allclose(h,h.conj().T,atol=1e-12)
        np.testing.assert_allclose(np.linalg.eigvalsh(h),
                                   [-xi,-xi,xi/2,xi/2,xi/2,xi/2],atol=1e-10)

    def test_majority_spin_lower_energy(self):
        S=np.array([[1.,.3],[.3,1.]],complex)
        H0=-5*S
        shifts=np.array([.5,.5])
        H,B=spinor_problem(H0,S,shifts,np.zeros((4,4)))
        energies=np.linalg.eigvalsh(np.linalg.solve(B,H))
        np.testing.assert_allclose(np.sort(energies),[-5.5,-5.5,-4.5,-4.5],
                                   atol=1e-11)
        np.testing.assert_allclose(H,H.conj().T,atol=1e-11)
        np.testing.assert_allclose(B,B.conj().T,atol=1e-11)

    def test_fermi_electron_conservation(self):
        energies=np.array([[-1,-1,1,1],[-1,-1,1,1.]])
        weights=np.array([.25,.75])
        occ,mu=fermi_dirac(energies,weights,2.0,.05)
        self.assertLess(abs(np.sum(occ*weights[:,None])-2.0),1e-11)
        np.testing.assert_allclose(occ[0],occ[1],atol=1e-11)

    def test_invalid_options(self):
        j=read_spin_job(ROOT/"examples"/"Fe_wire_spin.spin")
        with self.assertRaisesRegex(ValueError,"Choose spin"):
            replace(j.settings,spin=False,soc=False).validate()
        with self.assertRaisesRegex(ValueError,"mixing"):
            replace(j.settings,mixing=1.3).validate()
        with self.assertRaisesRegex(ValueError,"SOC requires"):
            replace(j.settings,spin=False,soc=True,xi={}).validate()


class TestFeWire(unittest.TestCase):
    def test_stoner_self_consistent_moment(self):
        j,r,ref,bands,polar=case("spin")
        self.assertTrue(r["converged"])
        self.assertLess(r["final_residual"],2e-5)
        self.assertGreater(r["total_moment_muB"],1.5)
        self.assertLess(r["total_moment_muB"],3.0)
        self.assertAlmostEqual(sum(r["atomic_charges"]),0,places=8)
        self.assertAlmostEqual(r["electrons"],
            np.sum(r["occupations"]*r["k_weights"][:,None]),places=9)
        self.assertEqual(bands.shape,(76,18))

    def test_pure_spin_channels_without_soc(self):
        _,_,_,bands,polar=case("spin")
        self.assertTrue(np.allclose(abs(polar),1,atol=1e-9))
        self.assertTrue(np.any(polar>0.99))
        self.assertTrue(np.any(polar<-.99))

    def test_soc_nonmagnetic_kramers(self):
        _,r,_,bands,polar=case("soc")
        self.assertTrue(r["converged"])
        self.assertEqual(r["iterations"],1)
        self.assertLess(abs(r["total_moment_muB"]),1e-8)
        # Wire uniatomica centrosimetrica: P+T garante duplas em TODO k.
        self.assertLess(np.max(abs(bands[:,0::2]-bands[:,1::2])),1e-8)

    def test_soc_stoner_mixed_spin(self):
        _,r,_,bands,polar=case("both")
        self.assertTrue(r["converged"])
        self.assertGreater(r["total_moment_muB"],1.5)
        self.assertLess(abs(r["total_moment_muB"]-
                            case("spin")[1]["total_moment_muB"]),.03)
        self.assertGreater(np.sum(abs(polar)<.8),0)
        self.assertGreater(np.max(abs(bands-case("spin")[3])),.005)
        self.assertLess(r["final_residual"],2e-5)

    def test_molecule_xyz_spin_soc(self):
        job=read_spin_job(ROOT/"examples"/"C_atom_both.spin")
        result,ref=scf(job.structure,job.kpoints,job.weights,job.settings,
                       quadrature=job.quadrature,cutoff=job.cutoff)
        self.assertTrue(result["converged"])
        self.assertGreater(result["total_moment_muB"],1.8)
        self.assertEqual(result["n_orbitals"],4)
        self.assertAlmostEqual(sum(result["atomic_charges"]),0.,places=9)
        self.assertEqual(len(result["integration_energies_eV"][0]),8)

    def test_isolated_atom_soc_without_stoner(self):
        job=read_spin_job(ROOT/"examples"/"C_atom_both.spin")
        settings=replace(job.settings,spin=False,soc=True)
        result,_=scf(job.structure,job.kpoints,job.weights,settings,
                     quadrature=job.quadrature,cutoff=job.cutoff)
        self.assertTrue(result["converged"])
        self.assertLess(abs(result["total_moment_muB"]),1e-8)
        self.assertLess(max(abs(result["atomic_charges"])),1e-9)

    def test_explicit_scf_nonconvergence(self):
        job=read_spin_job(ROOT/"examples"/"Fe_wire_spin.spin")
        settings=replace(job.settings,max_iter=1,workers=1)
        result,_=scf(job.structure,job.kpoints,job.weights,settings,
                     quadrature=job.quadrature,cutoff=job.cutoff)
        self.assertFalse(result["converged"])
        self.assertGreater(result["final_residual"],settings.tolerance)

    def test_spin_dos_normalization(self):
        job,res,ref,bands,polar=case("both")
        grid,up,down=spin_dos(res["integration_energies_eV"],
                              res["sz_integration"],job.weights,
                              sigma=.13,points=900)
        from scipy.integrate import trapezoid
        total=trapezoid(up+down,grid)
        self.assertAlmostEqual(total,18,places=4)
        self.assertTrue(np.all(up>=0))
        self.assertTrue(np.all(down>=0))

    def test_spin_flip_no_soc(self):
        job,res,ref,bands,polar=case("spin")
        settings=replace(job.settings,initial_moment=-2.5,workers=1)
        flipped,_=scf(job.structure,job.kpoints,job.weights,settings,
                      quadrature=job.quadrature,cutoff=job.cutoff)
        self.assertTrue(flipped["converged"])
        self.assertAlmostEqual(flipped["total_moment_muB"],
                               -res["total_moment_muB"],places=4)
        self.assertAlmostEqual(flipped["stoner_model_energy_eV"],
                               res["stoner_model_energy_eV"],places=7)

    def test_parallelism_same_result(self):
        job,res,ref,bands,polar=case("spin")
        settings=replace(job.settings,workers=1)
        serial,_=scf(job.structure,job.kpoints,job.weights,settings,
                     quadrature=job.quadrature,cutoff=job.cutoff)
        self.assertTrue(serial["converged"])
        self.assertAlmostEqual(serial["total_moment_muB"],res["total_moment_muB"],
                               places=7)
        self.assertAlmostEqual(serial["stoner_model_energy_eV"],
                               res["stoner_model_energy_eV"],places=7)

    def test_zero_soc_limit(self):
        job,r0,ref,bands,polar=case("spin")
        settings=replace(job.settings,soc=True,xi={"FE":{"d":1.e-9}})
        near_zero,ref2=scf(job.structure,job.kpoints,job.weights,settings,
                           quadrature=job.quadrature,cutoff=job.cutoff)
        bands0,pol=evaluate_bands(job.structure,job.band_path,ref2)
        self.assertTrue(near_zero["converged"])
        np.testing.assert_allclose(bands0,bands,atol=1e-7)

    def test_png_comparison_generated(self):
        cases={}
        for name in ("spin","soc","both"):
            job,res,ref,bands,polar=case(name)
            cases[name]={"bands_eV":bands,"spin_z_bands":polar,
                         "fermi_eV":res["fermi_eV"],
                         "special_kpoints":job.special_labels}
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/"comparacao.png"
            comparison_plot(cases,"Autes et al. 2006",p)
            zoom=Path(temp)/"zoom.png"
            plot_soc_zoom(cases["spin"],cases["both"],zoom)
            self.assertTrue(zoom.is_file())
            self.assertGreater(zoom.stat().st_size,10000)
            self.assertEqual(p.read_bytes()[:8],b"\x89PNG\r\n\x1a\n")


if __name__=="__main__":
    unittest.main()
