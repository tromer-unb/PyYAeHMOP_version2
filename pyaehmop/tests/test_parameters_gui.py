"""Regression tests for custom YAeHMOP parameter files and the desktop GUI."""
from dataclasses import replace
from pathlib import Path
import json
import os
import subprocess
import sys
import tempfile
import unittest
import numpy as np

from pyaehmop.basis import read_parameter_table
from pyaehmop.friendly import read_easy_job,calculate_easy_job
from pyaehmop.parameter_files import parse_parameter_file,load_parameter_selection
from pyaehmop.gui_logic import GuiSettings,build_config,load_config,render_config_to_file
from pyaehmop.gui_preview import load_structure_preview
from pyaehmop.gui_visuals import (
    draw_structure,draw_bands,draw_dos,draw_charges,draw_scf,format_summary)
from pyaehmop.gui_runner import run_gui_calculation

ROOT=Path(__file__).resolve().parents[1]


class TestParameterFiles(unittest.TestCase):
    def test_default_is_no_override(self):
        self.assertEqual(load_parameter_selection("default",ROOT),{})
        job=read_easy_job(ROOT/"examples_friendly/water_default_parameters.job")
        original=read_easy_job(ROOT/"examples_friendly/water_full.job")
        a=calculate_easy_job(job)
        b=calculate_easy_job(original)
        self.assertAlmostEqual(a["electronic_energy_eV"],
                               b["electronic_energy_eV"],places=11)
        np.testing.assert_allclose(a["mulliken_charges"],b["mulliken_charges"],atol=1e-10)

    def test_classic_parameters_file_same_default_result(self):
        job=read_easy_job(ROOT/"examples_friendly/water_custom_parameters.job")
        default=read_easy_job(ROOT/"examples_friendly/water_default_parameters.job")
        a=calculate_easy_job(job);b=calculate_easy_job(default)
        self.assertAlmostEqual(a["electronic_energy_eV"],
                               b["electronic_energy_eV"],places=11)
        np.testing.assert_allclose(a["mulliken_charges"],b["mulliken_charges"],atol=1e-10)

    def test_full_native_eht_parms_table(self):
        parsed=parse_parameter_file(ROOT/"eht_parms.dat")
        builtin=read_parameter_table()
        self.assertGreaterEqual(len(parsed),70)
        for element in ("H","O","FE","C","F"): 
            a=parsed[element];b=builtin[element]
            self.assertEqual(a["valence"],b["valence"])
            self.assertEqual(a["shells"],b["shells"])
        self.assertEqual(len(parsed["FE"]["shells"][-1].zetas),2)

    def test_full_table_equivalence_water(self):
        with tempfile.TemporaryDirectory() as tmp:
            src=ROOT/"examples_friendly/water.xyz"
            file=Path(tmp)/"calc.job"
            file.write_text(f"xyz\n{src}\nparameters={ROOT/'eht_parms.dat'}\n"
                            "calculations=energy,charges\n")
            a=calculate_easy_job(read_easy_job(file))
            b=calculate_easy_job(
                read_easy_job(ROOT/"examples_friendly/water_default_parameters.job"))
            self.assertAlmostEqual(a["electronic_energy_eV"],
                                   b["electronic_energy_eV"],places=11)

    def test_modified_parameter_changes_energy(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp);custom=path/"params.txt"
            custom.write_text("PARAMETERS\nH 1 1 1 1.3 -13.6\n"
                "O 8 6 2 2.275 -32.3 2 2.275 -15.5\nEND\n")
            job=path/"water.job"
            job.write_text(f"xyz\n{ROOT/'examples_friendly/water.xyz'}\n"
                           "parameters=params.txt\ncalculations=energy\n")
            report=calculate_easy_job(read_easy_job(job))
            ref=calculate_easy_job(
                read_easy_job(ROOT/"examples_friendly/water_default_parameters.job"))
            self.assertGreater(abs(report["electronic_energy_eV"]
                                   -ref["electronic_energy_eV"]),0.01)

    def test_missing_file_and_invalid_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                load_parameter_selection("missing_parameters.txt",tmp)
            path=Path(tmp)/"bad.txt"
            path.write_text("H 1 1 1 bad -13.6\n")
            with self.assertRaisesRegex(ValueError,"bad.txt"):
                parse_parameter_file(path)

    def test_custom_species_with_no_default_shells(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/"Ar.xyz").write_text("1\nArgon\nAr 0 0 0\n")
            (d/"Ar.par").write_text("AR 18 8 3 2.10 -20.0 3 1.80 -10.0\n")
            (d/"Ar.job").write_text("xyz\nAr.xyz\nparameters=Ar.par\n"
                                    "calculations=energy,charges\n")
            job=read_easy_job(d/"Ar.job")
            self.assertEqual(len(job.structure.atoms[0].shells),2)
            report=calculate_easy_job(job)
            self.assertEqual(report["orbitals"],4)
            self.assertEqual(report["electrons"],8)

    def test_spin_native_external_parameter_file(self):
        from pyaehmop.SPIN_TEST.config import read_spin_job
        ref=read_spin_job(ROOT/"SPIN_TEST/examples/Fe_wire_both.spin")
        ext=read_spin_job(
            ROOT/"SPIN_TEST/examples/Fe_wire_external_parameters.spin")
        self.assertEqual(ref.structure.atoms[0].shells,
                         ext.structure.atoms[0].shells)
        self.assertEqual(ref.structure.electrons,ext.structure.electrons)
        self.assertEqual(ref.settings,ext.settings)

    def test_unsupported_species_in_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)
            (path/"onlyH.txt").write_text("H 1 1 1 1.30 -13.60\n")
            (path/"O.xyz").write_text("1\nO\nO 0 0 0\n")
            (path/"O.job").write_text("xyz\nO.xyz\nparameters=onlyH.txt\n")
            with self.assertRaisesRegex(ValueError,"species"):
                read_easy_job(path/"O.job")


class TestGUIModel(unittest.TestCase):
    def test_default_job_is_native_friendly_file(self):
        source=ROOT/"examples_friendly/water.xyz"
        s=GuiSettings(structure_file=str(source),calculations=("energy","charges"))
        text=build_config(s)
        self.assertIn("parameters = default",text)
        self.assertIn("calculations = energy, charges",text)
        with tempfile.TemporaryDirectory() as tmp:
            file=render_config_to_file(s,tmp)
            loaded=read_easy_job(file)
            self.assertEqual(loaded.source_type,"xyz")
            self.assertEqual(set(loaded.tasks),{"energy","charges"})
            result=calculate_easy_job(loaded)
            self.assertAlmostEqual(result["electronic_energy_eV"],
                                   -162.5417635492,places=6)

    def test_legacy_inline_params_roundtrip(self):
        original=ROOT/"examples_friendly/diamond_cif.job"
        state=load_config(original)
        self.assertEqual(state.parameter_source,"inline")
        self.assertIn("parameters = C 6 4",build_config(state))

    def test_tk_window_construction(self):
        try:
            import tkinter as tk
            from pyaehmop.gui import PyYAeHMOPStudio
            root=tk.Tk()
        except (ImportError,tk.TclError if "tk" in locals() else ImportError) as exc:
            self.skipTest("No desktop Tk display available: "+str(exc))
        root.withdraw()
        try:
            ui=PyYAeHMOPStudio(root)
            root.update_idletasks()
            self.assertEqual(len(ui.tabs),7)
            self.assertIn("Bands",ui.tabs)
            self.assertIn("Structure",ui.tabs)
        finally:
            root.destroy()

    def test_custom_params_roundtrip(self):
        s=GuiSettings(
            structure_file=str(ROOT/"examples_friendly/water.xyz"),
            parameter_source="file",
            parameter_file=str(ROOT/"examples_friendly/water_parameters.txt"),
            calculations=("energy","charges","dos"),
            workers="2",output_dir=str(ROOT/"results/test_gui"),
            plots=False)
        with tempfile.TemporaryDirectory() as tmp:
            path=render_config_to_file(s,tmp)
            out=load_config(path)
            self.assertEqual(out.parameter_source,"file")
            self.assertEqual(out.parameter_file,s.parameter_file)
            self.assertEqual(set(out.calculations),set(s.calculations))
            self.assertEqual(out.workers,"2")

    def test_spin_soc_config_and_roundtrip(self):
        source=ROOT/"SPIN_TEST/examples/Fe_wire.extxyz"
        state=GuiSettings(structure_format="ase",structure_file=str(source),
                          mode="both",
                          kmesh="24",band_path="G:0,0,0; X:0.5,0,0",
                          divisions="75",cutoff="5.6",
                          quadrature="20 18 12",images="3")
        text=build_config(state)
        self.assertIn("spin = yes",text)
        self.assertIn("soc = yes",text)
        with tempfile.TemporaryDirectory() as tmp:
            path=render_config_to_file(state,tmp)
            self.assertEqual(path.suffix,".spin")
            loaded=load_config(path)
            self.assertEqual(loaded.mode,"both")
            from pyaehmop.SPIN_TEST.config import read_spin_job
            job=read_spin_job(path)
            self.assertTrue(job.settings.spin)
            self.assertTrue(job.settings.soc)
            self.assertEqual(len(job.kpoints),24)

    def test_invalid_duplicate_advanced_setting(self):
        s=GuiSettings(structure_file=str(ROOT/"examples_friendly/water.xyz"),
                      extra_options="workers = 3")
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            build_config(s)

    def test_invalid_missing_parameters(self):
        s=GuiSettings(structure_file=str(ROOT/"examples_friendly/water.xyz"),
                      parameter_source="file",parameter_file="/no/such/file")
        with self.assertRaises(FileNotFoundError):
            build_config(s)

    def test_gui_runner_eht(self):
        with tempfile.TemporaryDirectory() as tmp:
            state=GuiSettings(structure_file=str(ROOT/"examples_friendly/water.xyz"),
                 output_dir=tmp,calculations=("energy","charges","dos"),
                 plots=False)
            path=render_config_to_file(state,tmp)
            report_path=run_gui_calculation(path)
            r=json.loads(report_path.read_text())
            self.assertIn("mulliken_charges",r)
            self.assertIn("dos",r)
            self.assertAlmostEqual(r["electronic_energy_eV"],-162.5417635492,places=7)

    def test_gui_runner_soc_atom(self):
        with tempfile.TemporaryDirectory() as tmp:
            state=GuiSettings(structure_format="xyz",
                 structure_file=str(ROOT/"SPIN_TEST/examples/C_atom.xyz"),
                 mode="both",output_dir=tmp,workers="1",
                 spin_shell_parameters="stoner_c_p = 1.2\nsoc_c_p = 0.04",
                 plots=False,quadrature="12 12 8")
            path=render_config_to_file(state,tmp)
            data=json.loads(run_gui_calculation(path).read_text())
            self.assertTrue(data["converged"])
            self.assertGreater(data["total_moment_muB"],1.8)

    def test_structure_preview_and_matplotlib_figures(self):
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.figure import Figure
        view=load_structure_preview(ROOT/"examples_friendly/water.xyz","xyz")
        self.assertEqual(len(view.positions),3)
        fig=Figure(figsize=(6,4))
        draw_structure(fig,view)
        fig.canvas.draw()
        job=read_easy_job(ROOT/"examples_friendly/water_full.job")
        result=calculate_easy_job(job)
        for function in (draw_bands,draw_dos,draw_charges,draw_scf):
            function(fig,result)
            fig.canvas.draw()
        self.assertIn("Electronic energy",format_summary(result))


class TestAtomicParameterEditor(unittest.TestCase):
    def test_editor_exports_fe_double_zeta(self):
        from pyaehmop.parameter_editor_model import (
            default_rows,save_native_table,file_rows)
        rows=default_rows(["FE"])
        self.assertEqual(len(rows),3)
        self.assertEqual(next(r for r in rows if r["shell"]=="d")["n_zeta"],"2")
        with tempfile.TemporaryDirectory() as tmp:
            path=save_native_table(rows,Path(tmp)/"fe_native.txt")
            loaded=file_rows(path)
            self.assertEqual(rows,loaded)
            self.assertEqual(len(parse_parameter_file(path)["FE"]["shells"][-1].zetas),2)

    def test_editor_changes_onsite_energy(self):
        from pyaehmop.parameter_editor_model import (
            default_rows,save_native_table)
        rows=default_rows(["H","O"])
        oxygen=next(r for r in rows if r["symbol"]=="O" and r["shell"]=="p")
        oxygen["energy_eV"]="-16.1"
        with tempfile.TemporaryDirectory() as tmp:
            path=save_native_table(rows,Path(tmp)/"custom_native.txt")
            config=Path(tmp)/"water.job"
            config.write_text("xyz\n"+str(ROOT/"examples_friendly/water.xyz")+
                              "\nparameters="+str(path)+"\ncalculations=energy\n")
            new_energy=calculate_easy_job(read_easy_job(config))["electronic_energy_eV"]
            default_energy=calculate_easy_job(read_easy_job(
                ROOT/"examples_friendly/water_default_parameters.job"))["electronic_energy_eV"]
            self.assertGreater(abs(new_energy-default_energy),0.01)

    def test_editor_rejects_duplicates_and_inactive_shells(self):
        from pyaehmop.parameter_editor_model import (
            default_rows,native_table_text)
        rows=default_rows(["H"])
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            native_table_text(rows+rows)
        rows[0]["zeta1"]="0.0"
        with self.assertRaisesRegex(ValueError,"disabled"):
            native_table_text(rows)

    def test_editor_window_and_add_shell(self):
        try:
            import tkinter as tk
            from pyaehmop.parameter_editor_dialog import AtomicParameterEditor
            root=tk.Tk()
        except (ImportError,tk.TclError if "tk" in locals() else ImportError) as exc:
            self.skipTest("No Tk display: "+str(exc))
        root.withdraw()
        try:
            dlg=AtomicParameterEditor(root,["C"])
            self.assertEqual(len(dlg.rows),2)
            dlg._add_row()
            self.assertEqual(len(dlg.rows),3)
            self.assertEqual(dlg.rows[-1]["shell"],"d")
            dlg.vars["energy_eV"].set("-4.5")
            dlg._update_row()
            self.assertEqual(dlg.rows[-1]["energy_eV"],"-4.5")
            dlg.top.destroy()
        finally:
            root.destroy()


if __name__=="__main__":
    unittest.main()
