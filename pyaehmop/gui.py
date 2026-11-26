"""PyYAeHMOP Studio — native desktop GUI for friendly EHT and spin/SOC.

Launch: python3 -m pyaehmop.gui  (or pyaehmop-gui after installation)
Requires a graphical session, tkinter, matplotlib, and ASE for non-YAeHMOP
structure previews. The solver always runs in a separate Python process.
"""
import json
import os
from pathlib import Path
import queue
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from dataclasses import fields

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from .gui_logic import GuiSettings, TASKS, build_config, load_config, render_config_to_file
from .gui_preview import load_structure_preview
from .gui_visuals import (draw_structure,draw_bands,draw_dos,draw_charges,
                          draw_scf,format_summary)


MODE_LABELS={
    "Standard EHT":"eht",
    "Spin SCF (↑ / ↓)":"spin",
    "Spin-orbit only (SOC)":"soc",
    "Spin SCF + SOC":"both",
}
MODE_NAMES={code:label for label,code in MODE_LABELS.items()}

FILE_TYPES=[
    ("Structure / input files","*.xyz *.extxyz *.cif *.vasp *.POSCAR *.CONTCAR *.in *.input *.txt"),
    ("All files","*.*")
]


class PyYAeHMOPStudio:
    def __init__(self,root):
        self.root=root
        self.root.title("PyYAeHMOP Studio  •  EHT / Spin / SOC")
        self.root.geometry("1440x860")
        self.root.minsize(1030,685)
        self.root.configure(background="#f1f4f8")
        self.messages=queue.Queue()
        self.process=None
        self.busy=False
        self.result=None
        self.view=None
        self.figures={}
        self.canvases={}
        self.tabs={}
        self.settings=GuiSettings()
        self._configure_style()
        self._layout()
        self._apply_state(self.settings)
        self.root.protocol("WM_DELETE_WINDOW",self._close)
        self.root.after(120,self._poll_messages)

    def _configure_style(self):
        style=ttk.Style(self.root)
        try:style.theme_use("clam")
        except tk.TclError:pass
        style.configure("TFrame",background="#f5f7fb")
        style.configure("TLabel",background="#f5f7fb",foreground="#283649",
                        font=("DejaVu Sans",9))
        style.configure("TLabelframe",background="#f5f7fb",
                        bordercolor="#d6deea",relief="flat")
        style.configure("TLabelframe.Label",background="#f5f7fb",
                        foreground="#214b79",font=("DejaVu Sans",10,"bold"))
        style.configure("TButton",font=("DejaVu Sans",9),padding=(9,7))
        style.configure("Accent.TButton",background="#116b9d",foreground="white",
                        bordercolor="#116b9d",font=("DejaVu Sans",10,"bold"),padding=(15,8))
        style.map("Accent.TButton",background=[("active","#0b5883")])
        style.configure("TEntry",padding=5)
        style.configure("TCombobox",padding=4)
        style.configure("TCheckbutton",background="#f5f7fb",foreground="#243447")
        style.configure("TNotebook",background="#e9eff7")
        style.configure("TNotebook.Tab",padding=(15,8),font=("DejaVu Sans",9))
        style.map("TNotebook.Tab",background=[("selected","#ffffff")])
        style.configure("Horizontal.TProgressbar",background="#167caa")

    def _layout(self):
        header=tk.Frame(self.root,bg="#132940",height=69)
        header.pack(side="top",fill="x")
        header.pack_propagate(False)
        tk.Label(header,text="◈  PyYAeHMOP STUDIO",bg="#132940",fg="#f5fbff",
                 font=("DejaVu Sans",16,"bold")).pack(side="left",padx=(24,18),pady=15)
        tk.Label(header,text="Extended Hückel Theory   •   Magnetism   •   Spin–Orbit",
                 bg="#132940",fg="#a8cbdf",
                 font=("DejaVu Sans",10)).pack(side="left",pady=15)
        self.run_button=ttk.Button(header,text="▶  Run calculation",style="Accent.TButton",
                                   command=self._run)
        self.run_button.pack(side="right",padx=(10,22),pady=12)
        self.stop_button=ttk.Button(header,text="Stop",command=self._stop,state="disabled")
        self.stop_button.pack(side="right",pady=12)
        actions=ttk.Frame(self.root,padding=(13,8))
        actions.pack(fill="x")
        for text,callback in (
            ("Open structure",self._browse_structure),
            ("Open .job / .spin",self._open_input),
            ("Save input",self._save_input),
            ("Preview structure",self._preview),
            ("Open results JSON",self._open_results),
        ):
            ttk.Button(actions,text=text,command=callback).pack(side="left",padx=4)
        self.status=tk.StringVar(value="Ready — select a structure file to begin")
        tk.Label(actions,textvariable=self.status,bg="#f5f7fb",fg="#486581",
                 font=("DejaVu Sans",9)).pack(side="right",padx=12)
        frame=ttk.Frame(self.root)
        frame.pack(fill="both",expand=True,padx=10,pady=(0,6))
        paned=ttk.Panedwindow(frame,orient="horizontal")
        paned.pack(fill="both",expand=True)
        left=ttk.Frame(paned,width=405)
        right=ttk.Frame(paned)
        paned.add(left,weight=0)
        paned.add(right,weight=1)
        self._setup_scroll(left)
        self._setup_views(right)
        foot=ttk.Frame(self.root,padding=(16,4))
        foot.pack(side="bottom",fill="x")
        self.progress=ttk.Progressbar(foot,mode="indeterminate",length=145)
        self.progress.pack(side="right")
        ttk.Label(foot,text="Computations are performed by the existing Python EHT engine.",
                  foreground="#68778b").pack(side="left")

    def _setup_scroll(self,host):
        canvas=tk.Canvas(host,bg="#f5f7fb",highlightthickness=0,width=400)
        scrollbar=ttk.Scrollbar(host,orient="vertical",command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left",fill="both",expand=True)
        scrollbar.pack(side="right",fill="y")
        content=ttk.Frame(canvas,padding=(9,10))
        window=canvas.create_window((0,0),window=content,anchor="nw")
        content.bind("<Configure>",lambda e:canvas.configure(
            scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>",lambda e:canvas.itemconfigure(window,width=e.width))
        canvas.bind("<Enter>",lambda e:canvas.bind_all(
            "<MouseWheel>",lambda ev:canvas.yview_scroll(
                -1*(ev.delta//120),"units")))
        canvas.bind("<Leave>",lambda e:canvas.unbind_all("<MouseWheel>"))
        self.controls=content
        self.vars={}
        for f in fields(GuiSettings):
            if f.name not in ("calculations","mode","spin_shell_parameters","extra_options"):
                self.vars[f.name]=(tk.BooleanVar() if f.name=="plots" else tk.StringVar())
        self.field_widgets={}
        self.task_vars={key:tk.BooleanVar(value=False) for key in ("all",*TASKS)}
        self.mode_var=tk.StringVar(value=MODE_NAMES["eht"])
        self._input_controls()
        self._task_controls()
        self._scientific_controls()
        self._spin_controls()
        self._output_controls()
        self._advanced_controls()

    def _group(self,title):
        group=ttk.LabelFrame(self.controls,text=title,padding=(10,8))
        group.pack(fill="x",padx=(3,7),pady=(0,11))
        return group

    def _field(self,host,label,key,browse=None,values=None):
        row=ttk.Frame(host)
        row.pack(fill="x",pady=3)
        ttk.Label(row,text=label,width=16).pack(side="left")
        if values is not None:
            control=ttk.Combobox(row,textvariable=self.vars[key],values=values,
                                 state="readonly",width=17)
        else:
            control=ttk.Entry(row,textvariable=self.vars[key],width=24)
        control.pack(side="left",fill="x",expand=True,padx=(4,3))
        self.field_widgets[key]=control
        if browse is not None:
            ttk.Button(row,text="…",width=3,command=browse).pack(side="right")
        return control

    def _input_controls(self):
        group=self._group("01   Structure & parameters")
        self._field(group,"Structure format","structure_format",
                    values=("yaehmop","xyz","cif","ase"))
        self._field(group,"Structure file","structure_file",self._browse_structure)
        ttk.Label(group,text="EHT basis parameters").pack(anchor="w",pady=(8,0))
        self.parameter_mode=ttk.Combobox(
            group,textvariable=self.vars["parameter_source"],
            values=("default","file","inline"),state="readonly")
        self.parameter_mode.pack(fill="x",pady=(2,4))
        self.parameter_mode.bind("<<ComboboxSelected>>",lambda ev:self._sync_parameters())
        self.parameter_entry=self._field(group,"Parameter file / line","parameter_file",
                                         self._browse_parameters)
        ttk.Button(group,text="Edit atomic shell table…",
                   command=self._edit_parameters).pack(anchor="w",pady=(4,4))
        ttk.Label(group,text="Supports native YAeHMOP eht_parms.dat and classic PARAMETERS.",
                  foreground="#56687a",wraplength=338).pack(anchor="w",pady=(4,0))
        ttk.Label(group,text="Calculation mode").pack(anchor="w",pady=(10,2))
        self.mode_input=ttk.Combobox(group,textvariable=self.mode_var,
                                     values=tuple(MODE_LABELS),state="readonly")
        self.mode_input.pack(fill="x")
        self.mode_input.bind("<<ComboboxSelected>>",lambda ev:self._sync_mode())

    def _task_controls(self):
        self.task_frame=self._group("02   Requested calculations")
        self.task_widgets=[]
        options=[("All available","all"),("Electronic energy","energy"),
                 ("Mulliken charges","charges"),("Band structure","bands"),
                 ("Density of states","dos"),("COOP / COHP bonds","bonds"),
                 ("Fragment orbitals","fmo"),("Orbital levels","orbitals"),
                 ("Matrices","matrices")]
        for i,(label,key) in enumerate(options):
            item=ttk.Checkbutton(self.task_frame,text=label,
                    variable=self.task_vars[key],command=lambda k=key:self._sync_all(k))
            item.grid(row=i//2,column=i%2,sticky="w",pady=4,padx=2)
            self.task_widgets.append(item)
        self.task_note=ttk.Label(self.task_frame,
            text="For spin/SOC, the existing solver always computes spin DOS; "
                 "bands are generated when a band path is supplied.",
            wraplength=330,foreground="#68778b")
        self.task_note.grid(row=5,column=0,columnspan=2,sticky="w",pady=(5,2))

    def _scientific_controls(self):
        group=self._group("03   Engine, k points & DOS")
        for label,key in (
            ("CPU workers","workers"),
            ("STO quadrature","quadrature"),
            ("Cutoff (Å)","cutoff"),
            ("Periodic images","images"),
            ("K-point grid","kmesh"),
            ("Explicit k points","kpoints"),
            ("K-grid shift","k_shift"),
            ("Band path","band_path"),
            ("Band divisions","divisions"),
            ("DOS width (eV)","dos_sigma"),
            ("DOS points","dos_points"),
            ("Charge","charge"),
            ("Electron count","electrons"),
            ("Coupling K","k"),
            ("Weighted EHT","weighted"),
            ("Bond pairs","bonds"),
            ("FMO fragments","fragments"),
            ("ASE file format","ase_format"),
            ("Trajectory frame","frame"),
        ):
            self._field(group,label,key)
        ttk.Label(group,text="Leave kmesh blank for default integration; "
                  "use band_path = auto for ASE high-symmetry paths.",
                  wraplength=325,foreground="#68778b").pack(anchor="w",pady=(5,0))

    def _spin_controls(self):
        self.spin_frame=self._group("04   Spin / spin–orbit settings")
        for label,key in (
            ("Initial moment","initial_moment"),
            ("SCF mixing","mixing"),
            ("kBT (eV)","temperature_ev"),
            ("SCF tolerance","tolerance"),
            ("Max iterations","max_iter"),
        ):
            self._field(self.spin_frame,label,key)
        ttk.Label(self.spin_frame,text="Onsite Stoner & SOC shell strengths (eV)").pack(
            anchor="w",pady=(7,3))
        self.spin_text=tk.Text(self.spin_frame,height=3,font=("DejaVu Sans Mono",9),
                               wrap="none",borderwidth=1,relief="solid")
        self.spin_text.pack(fill="x")
        ttk.Label(self.spin_frame,
            text="Example: stoner_fe_d = 0.88 and soc_fe_d = 0.06. "
                 "These are demonstration parameters, not a universal fit.",
            wraplength=335,foreground="#65748a").pack(anchor="w",pady=(4,0))

    def _output_controls(self):
        group=self._group("05   Output settings")
        self._field(group,"Output directory","output_dir",self._browse_output)
        self._field(group,"File formats","outputs")
        ttk.Checkbutton(group,text="Generate PNG plots",
                        variable=self.vars["plots"]).pack(anchor="w",pady=(5,1))
        ttk.Label(group,text="Results open automatically in the tabs on the right.",
                  foreground="#65748a",wraplength=325).pack(anchor="w",pady=(2,0))

    def _advanced_controls(self):
        group=self._group("06   Advanced input options")
        ttk.Label(group,text="Additional key = value lines (optional)").pack(anchor="w")
        self.advanced_text=tk.Text(group,height=3,font=("DejaVu Sans Mono",9),
                                   wrap="none",borderwidth=1,relief="solid")
        self.advanced_text.pack(fill="x",pady=4)
        ttk.Label(group,text="Unrecognized options are validated by the engine.",
                  wraplength=320,foreground="#65748a").pack(anchor="w",pady=(0,2))

    def _setup_views(self,parent):
        self.notebook=ttk.Notebook(parent)
        self.notebook.pack(fill="both",expand=True,padx=(2,4),pady=(0,1))
        for label in ("Structure","Bands","DOS","Charges","SCF"):
            tab=ttk.Frame(self.notebook)
            self.notebook.add(tab,text=label)
            self.tabs[label]=tab
            fig=Figure(figsize=(8.5,5.6),dpi=100)
            self.figures[label]=fig
            canvas=FigureCanvasTkAgg(fig,master=tab)
            widget=canvas.get_tk_widget()
            widget.pack(fill="both",expand=True)
            self.canvases[label]=canvas
            toolbar=NavigationToolbar2Tk(canvas,tab,pack_toolbar=False)
            toolbar.update()
            toolbar.pack(side="bottom",fill="x")
            ax=fig.add_subplot(111)
            ax.set_axis_off()
            ax.text(.5,.5,"Open a structure or run a calculation to display results.",
                    transform=ax.transAxes,ha="center",va="center",fontsize=12,
                    color="#688098")
            canvas.draw_idle()
        summary=ttk.Frame(self.notebook)
        self.notebook.add(summary,text="Summary")
        self.tabs["Summary"]=summary
        self.summary_text=tk.Text(summary,wrap="word",font=("DejaVu Sans Mono",10),
                                  bg="#ffffff",fg="#27384d",padx=18,pady=18)
        self.summary_text.pack(fill="both",expand=True)
        self.summary_text.insert("1.0","Select a structure, configure the model and click Run.")
        console=ttk.Frame(self.notebook)
        self.notebook.add(console,text="Console")
        self.tabs["Console"]=console
        self.console=tk.Text(console,wrap="word",font=("DejaVu Sans Mono",9),
                             bg="#111d2c",fg="#e3f1f8",insertbackground="#ffffff",
                             padx=15,pady=14)
        self.console.pack(fill="both",expand=True)
        self._append_log("PyYAeHMOP Studio started. Original scientific engine retained.")

    def _apply_state(self,state):
        self.settings=state
        for key,var in self.vars.items():
            if hasattr(state,key):
                var.set(getattr(state,key))
        for key,var in self.task_vars.items():
            var.set(key in state.calculations)
        self.mode_var.set(MODE_NAMES.get(state.mode,MODE_NAMES["eht"]))
        self.spin_text.configure(state="normal")
        self.spin_text.delete("1.0","end")
        self.spin_text.insert("1.0",state.spin_shell_parameters)
        self.advanced_text.delete("1.0","end")
        self.advanced_text.insert("1.0",state.extra_options)
        self._sync_all("all",notify=False)
        self._sync_parameters()
        self._sync_mode()

    def _collect(self):
        data={key:var.get() for key,var in self.vars.items()}
        data["mode"]=MODE_LABELS[self.mode_var.get()]
        selected=tuple(k for k,v in self.task_vars.items() if v.get())
        data["calculations"]=selected
        data["spin_shell_parameters"]=self.spin_text.get("1.0","end").strip()
        data["extra_options"]=self.advanced_text.get("1.0","end").strip()
        return GuiSettings(**data)

    def _sync_all(self,key,notify=True):
        if key=="all" and self.task_vars["all"].get():
            for other in TASKS:self.task_vars[other].set(False)
        elif key!="all" and self.task_vars[key].get():
            self.task_vars["all"].set(False)
        if notify:
            self.status.set("Calculation selection updated")

    def _sync_parameters(self):
        editable=self.vars["parameter_source"].get() in ("file","inline")
        self.parameter_entry.configure(state="normal" if editable else "disabled")

    def _sync_mode(self):
        mode=MODE_LABELS.get(self.mode_var.get(),"eht")
        for widget in self.task_widgets:
            widget.state(["!disabled"] if mode=="eht" else ["disabled"])
        self.task_note.configure(
            text=("Choose properties to compute for standard extended Hückel."
                  if mode=="eht" else
                  "Spin/SOC automatically computes spin-resolved DOS and SCF; "
                  "include a band path to produce spin-projected bands."))
        for item in self.spin_frame.winfo_children():
            # Avoid disabling labels to keep instructions readable.
            if isinstance(item,ttk.Frame):
                for child in item.winfo_children():
                    if isinstance(child,ttk.Entry):
                        child.configure(state="normal" if mode!="eht" else "disabled")
        self.spin_text.configure(state="normal" if mode!="eht" else "disabled")
        for key in ("k","weighted","bonds","fragments","kpoints","k_shift","frame"):
            self.field_widgets[key].configure(state="normal" if mode=="eht" else "disabled")

    def _browse_structure(self):
        file=filedialog.askopenfilename(parent=self.root,title="Choose structure file",
                                       filetypes=FILE_TYPES)
        if not file:return
        self.vars["structure_file"].set(file)
        ext=Path(file).suffix.lower()
        kind=("cif" if ext==".cif" else "xyz" if ext==".xyz"
              else "ase" if ext in (".extxyz",".vasp",".traj",".poscar")
              or Path(file).name.upper() in ("POSCAR","CONTCAR")
              else "yaehmop")
        self.vars["structure_format"].set(kind)
        self._preview()

    def _browse_parameters(self):
        file=filedialog.askopenfilename(parent=self.root,
            title="Select YAeHMOP parameter table",
            filetypes=[("Parameter files","*.txt *.dat *.par *.params"),("All files","*.*")])
        if file:
            self.vars["parameter_source"].set("file")
            self.vars["parameter_file"].set(file)
            self._sync_parameters()

    def _edit_parameters(self):
        """Open an interactive shell editor; never overwrite the defaults."""
        try:
            from .parameter_editor_dialog import show_atomic_parameter_editor
            state=self._collect()
            species=[]
            if state.structure_file and Path(state.structure_file).is_file():
                view=(self.view if self.view is not None else
                      load_structure_preview(state.structure_file,
                          state.structure_format,state.ase_format,state.frame or "0"))
                species=view.symbols
            param_file=(state.parameter_file if state.parameter_source=="file"
                        else None)
            def selected(path):
                self.vars["parameter_source"].set("file")
                self.vars["parameter_file"].set(path)
                self._sync_parameters()
                self._append_log("Selected custom YAeHMOP atomic table: "+path)
                self.status.set("Edited atomic parameter table selected")
            show_atomic_parameter_editor(self.root,species,param_file,selected)
        except Exception as exc:
            messagebox.showerror("Atomic parameter editor",str(exc),parent=self.root)

    def _browse_output(self):
        directory=filedialog.askdirectory(parent=self.root,title="Choose output directory")
        if directory:self.vars["output_dir"].set(directory)

    def _open_input(self):
        file=filedialog.askopenfilename(parent=self.root,
            title="Open a friendly .job or .spin file",
            filetypes=[("PyYAeHMOP inputs","*.job *.spin"),("All files","*.*")])
        if not file:return
        try:
            self._apply_state(load_config(file))
            self._append_log("Loaded input: "+file)
            self.status.set("Input file loaded")
            self._preview()
        except Exception as exc:
            messagebox.showerror("Could not open input",str(exc),parent=self.root)

    def _save_input(self):
        try:
            state=self._collect()
            default_ext=".job" if state.mode=="eht" else ".spin"
            text=build_config(state)
            file=filedialog.asksaveasfilename(parent=self.root,title="Save calculation input",
                defaultextension=default_ext,
                filetypes=[("PyYAeHMOP input","*"+default_ext),("All files","*.*")])
            if not file:return
            path=Path(file)
            path.write_text(text,encoding="utf-8")
            self._append_log("Saved input: "+str(path))
            self.status.set("Input saved")
        except Exception as exc:
            messagebox.showerror("Invalid input",str(exc),parent=self.root)

    def _preview(self):
        state=self._collect()
        if not state.structure_file.strip():
            self.status.set("Choose a structure file to preview")
            return
        try:
            self.view=load_structure_preview(state.structure_file,
                 state.structure_format,state.ase_format,state.frame or "0")
            draw_structure(self.figures["Structure"],self.view)
            self.canvases["Structure"].draw_idle()
            self.notebook.select(self.tabs["Structure"])
            self.status.set(f"Loaded {len(self.view.symbols)} atoms • "
                            f"{len(self.view.vectors)} periodic directions")
        except Exception as exc:
            self._append_log("Structure preview error: "+str(exc))
            messagebox.showerror("Structure preview",str(exc),parent=self.root)

    def _open_results(self):
        file=filedialog.askopenfilename(parent=self.root,title="Open calculation JSON",
                  filetypes=[("JSON files","*.json"),("All files","*.*")])
        if file:
            self._display_result(file)

    def _display_result(self,file):
        try:
            report=json.loads(Path(file).read_text(encoding="utf-8"))
            self.result=report
            for key,plot in (("Bands",draw_bands),("DOS",draw_dos),
                             ("Charges",draw_charges),("SCF",draw_scf)):
                plot(self.figures[key],report)
                self.canvases[key].draw_idle()
            self.summary_text.delete("1.0","end")
            self.summary_text.insert("1.0",format_summary(report))
            self.notebook.select(self.tabs["Summary"])
            self.status.set("Results loaded: "+Path(file).name)
            self._append_log("Loaded results: "+str(file))
        except Exception as exc:
            messagebox.showerror("Could not load results",str(exc),parent=self.root)

    def _set_busy(self,flag):
        self.busy=flag
        self.run_button.configure(state="disabled" if flag else "normal")
        self.stop_button.configure(state="normal" if flag else "disabled")
        if flag:self.progress.start(10)
        else:self.progress.stop()

    def _run(self):
        if self.busy:return
        try:
            state=self._collect()
            output=(Path(state.output_dir).expanduser().resolve()
                    if state.output_dir else Path.home()/"pyaehmop_gui_results")
            path=render_config_to_file(state,output)
        except Exception as exc:
            messagebox.showerror("Invalid calculation setup",str(exc),parent=self.root)
            return
        self._set_busy(True)
        self.status.set("Calculation running in a separate process…")
        self._append_log("\nStarting calculation: "+str(path))
        self.notebook.select(self.tabs["Console"])
        thread=threading.Thread(target=self._worker,args=(path,),daemon=True)
        thread.start()

    def _worker(self,config):
        output_file=None
        lines=[]
        cmd=[sys.executable,"-m","pyaehmop.gui_runner",str(config)]
        env=os.environ.copy()
        env["OPENBLAS_NUM_THREADS"]="1"
        env["OMP_NUM_THREADS"]="1"
        try:
            process=subprocess.Popen(cmd,stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,text=True,bufsize=1,env=env)
            self.process=process
            for line in process.stdout:
                line=line.rstrip("\n")
                lines.append(line)
                self.messages.put(("line",line))
                if line.startswith("RESULT_FILE="):
                    output_file=line.split("=",1)[1].strip()
            code=process.wait()
            self.messages.put(("done",code,output_file,"\n".join(lines[-12:])))
        except Exception as exc:
            self.messages.put(("done",1,None,str(exc)))
        finally:
            self.process=None

    def _poll_messages(self):
        try:
            while True:
                item=self.messages.get_nowait()
                if item[0]=="line":
                    self._append_log(item[1])
                elif item[0]=="done":
                    _,code,output_file,details=item
                    self._set_busy(False)
                    if code==0 and output_file:
                        self._display_result(output_file)
                        self.status.set("Calculation completed successfully")
                    else:
                        self.status.set("Calculation failed — inspect Console")
                        self.notebook.select(self.tabs["Console"])
                        messagebox.showerror("Calculation failed",
                            details[-2400:] or f"Exit code {code}",parent=self.root)
        except queue.Empty:
            pass
        self.root.after(120,self._poll_messages)

    def _stop(self):
        p=self.process
        if p is not None and p.poll() is None:
            try:p.terminate()
            except OSError:pass
            self._append_log("Stop requested.")
            self.status.set("Stopping calculation…")

    def _append_log(self,line):
        self.console.insert("end",str(line)+"\n")
        self.console.see("end")

    def _close(self):
        if self.busy:
            ok=messagebox.askyesno("Calculation running",
                "Stop the calculation and close the window?",parent=self.root)
            if not ok:return
            self._stop()
        self.root.destroy()


def main():
    root=tk.Tk()
    PyYAeHMOPStudio(root)
    root.mainloop()


if __name__=="__main__":
    main()
