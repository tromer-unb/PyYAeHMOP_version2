"""Tkinter editor for atomic Slater shell parameters.

All controls operate on a copy of shell records. Saving creates a new
native YAeHMOP table; the bundled source and scientific code stay unchanged.
"""
from pathlib import Path
import tkinter as tk
from tkinter import ttk,filedialog,messagebox

from .parameter_editor_model import (
    COLUMNS,default_rows,file_rows,validate_row,save_native_table)


LABELS={
    "symbol":"Element",
    "atomic_number":"Atomic Z",
    "valence":"Valence e⁻",
    "n_zeta":"Nζ",
    "n":"n",
    "shell":"Shell",
    "energy_eV":"Hii (eV)",
    "zeta1":"ζ₁",
    "zeta2":"ζ₂",
    "coeff1":"c₁",
    "coeff2":"c₂",
}
WIDTHS={
    "symbol":82,"atomic_number":85,"valence":95,"n_zeta":52,
    "n":48,"shell":62,"energy_eV":100,
    "zeta1":90,"zeta2":90,"coeff1":85,"coeff2":85
}


class AtomicParameterEditor:
    def __init__(self,parent,selected_species=None,source_file=None,on_save=None):
        self.parent=parent
        self.on_save=on_save
        self.species=list(dict.fromkeys(str(s).upper() for s in
                           (selected_species or [])))
        self.rows=default_rows(self.species or None)
        self.top=tk.Toplevel(parent)
        self.top.title("YAeHMOP atomic parameter editor — Slater shells")
        self.top.geometry("1160x650")
        self.top.minsize(945,565)
        self.top.transient(parent)
        self.top.configure(background="#f5f7fb")
        self.vars={key:tk.StringVar() for key in COLUMNS}
        self.current=None
        self._build()
        if source_file and Path(source_file).is_file():
            try:
                self.rows=file_rows(source_file)
            except Exception as exc:
                messagebox.showwarning("Parameter file",
                    f"Could not load the selected parameter file: {exc}",
                    parent=self.top)
        self._refresh()
        self.top.grab_set()

    def _build(self):
        heading=ttk.Frame(self.top,padding=(18,13))
        heading.pack(fill="x")
        ttk.Label(heading,text="ATOMIC SLATER SHELLS",
                  font=("DejaVu Sans",14,"bold"),
                  foreground="#16395b").pack(anchor="w")
        ttk.Label(heading,
            text="Edit Hii, orbital exponent(s), and double-zeta coefficients. "
                 "Changes are saved as a separate native YAeHMOP table.",
            foreground="#60758b",wraplength=1090).pack(anchor="w",pady=(3,0))

        toolbar=ttk.Frame(self.top,padding=(16,0))
        toolbar.pack(fill="x",pady=(0,8))
        ttk.Button(toolbar,text="Load defaults for structure",
                   command=self._load_defaults).pack(side="left",padx=(0,7))
        ttk.Button(toolbar,text="Import YAeHMOP file…",
                   command=self._import_file).pack(side="left",padx=(0,7))
        ttk.Label(toolbar,text=("Species: "+", ".join(self.species)
                  if self.species else "All bundled species"),
                  foreground="#41688c").pack(side="right")

        table_frame=ttk.Frame(self.top,padding=(13,2))
        table_frame.pack(fill="both",expand=True)
        self.tree=ttk.Treeview(table_frame,columns=COLUMNS,show="headings",
                               selectmode="browse",height=12)
        for key in COLUMNS:
            self.tree.heading(key,text=LABELS[key])
            self.tree.column(key,width=WIDTHS[key],minwidth=43,stretch=True,anchor="center")
        vscroll=ttk.Scrollbar(table_frame,orient="vertical",command=self.tree.yview)
        hscroll=ttk.Scrollbar(table_frame,orient="horizontal",command=self.tree.xview)
        self.tree.configure(yscrollcommand=vscroll.set,xscrollcommand=hscroll.set)
        self.tree.grid(row=0,column=0,sticky="nsew")
        vscroll.grid(row=0,column=1,sticky="ns")
        hscroll.grid(row=1,column=0,sticky="ew")
        table_frame.columnconfigure(0,weight=1)
        table_frame.rowconfigure(0,weight=1)
        self.tree.bind("<<TreeviewSelect>>",self._select_row)

        editor=ttk.LabelFrame(self.top,text="Selected orbital shell",padding=(12,10))
        editor.pack(fill="x",padx=17,pady=(10,6))
        for index,key in enumerate(COLUMNS):
            row=index//6
            column=index%6
            panel=ttk.Frame(editor)
            panel.grid(row=row,column=column,sticky="ew",padx=5,pady=(2,8))
            ttk.Label(panel,text=LABELS[key]).pack(anchor="w")
            if key=="shell":
                widget=ttk.Combobox(panel,textvariable=self.vars[key],
                                    values=("s","p","d","f"),state="readonly",width=10)
            elif key=="n_zeta":
                widget=ttk.Combobox(panel,textvariable=self.vars[key],
                                    values=("1","2"),state="readonly",width=10)
            else:
                widget=ttk.Entry(panel,textvariable=self.vars[key],width=13)
            widget.pack(fill="x")
            editor.columnconfigure(column,weight=1)

        buttons=ttk.Frame(self.top,padding=(15,5))
        buttons.pack(fill="x")
        ttk.Button(buttons,text="Apply changes to selected shell",
                   command=self._update_row).pack(side="left",padx=4)
        ttk.Button(buttons,text="Add new shell",
                   command=self._add_row).pack(side="left",padx=4)
        ttk.Button(buttons,text="Remove selected shell",
                   command=self._delete_row).pack(side="left",padx=4)
        ttk.Button(buttons,text="Cancel",
                   command=self.top.destroy).pack(side="right",padx=(5,0))
        ttk.Button(buttons,text="Save YAeHMOP table and use",
                   style="Accent.TButton",
                   command=self._save).pack(side="right",padx=5)

    def _refresh(self,index=0):
        for item in self.tree.get_children():
            self.tree.delete(item)
        for i,row in enumerate(self.rows):
            self.tree.insert("", "end",iid=str(i),
                             values=tuple(row.get(key,"") for key in COLUMNS))
        self.current=None
        if self.rows:
            index=min(max(index,0),len(self.rows)-1)
            self.tree.selection_set(str(index))
            self.tree.see(str(index))
            self._select_row()
        else:
            self._load_fields(self._blank_row())

    def _load_fields(self,row):
        for key in COLUMNS:
            self.vars[key].set(str(row.get(key,"")))

    def _blank_row(self):
        element=(self.species[0] if self.species else "C")
        try:
            from ase.data import atomic_numbers
            number=atomic_numbers.get(element.capitalize(),6)
        except ImportError:
            number=6
        return dict(zip(COLUMNS,(element,str(number),"4","1","2","s",
                                  "-10.0","1.5","0.0","1.0","0.0")))

    def _select_row(self,event=None):
        selected=self.tree.selection()
        if not selected:
            return
        index=int(selected[0])
        if not 0<=index<len(self.rows):
            return
        self.current=index
        self._load_fields(self.rows[index])

    def _read_fields(self):
        return validate_row({key:self.vars[key].get() for key in COLUMNS})

    def _update_row(self):
        if self.current is None:
            messagebox.showinfo("Select a shell",
                "Select a row first, or click Add new shell.",parent=self.top)
            return
        try:
            row=self._read_fields()
            self.rows[self.current]=row
            self._refresh(self.current)
        except Exception as exc:
            messagebox.showerror("Invalid shell",str(exc),parent=self.top)

    def _add_row(self):
        """Insert a new editable shell with safe initial values."""
        row=self._blank_row()
        if self.current is not None:
            previous=self.rows[self.current]
            for key in ("symbol","atomic_number","valence"):
                row[key]=previous[key]
        existing={r["shell"] for r in self.rows if r["symbol"]==row["symbol"]}
        for shell,quantum in (("s","2"),("p","2"),("d","3"),("f","4")):
            if shell not in existing:
                row["shell"]=shell
                row["n"]=quantum
                break
        self.rows.append(row)
        self._refresh(len(self.rows)-1)

    def _delete_row(self):
        if self.current is None:return
        if not messagebox.askyesno("Delete shell",
            "Remove the selected shell from the custom parameter table?",
            parent=self.top):
            return
        self.rows.pop(self.current)
        self._refresh()

    def _load_defaults(self):
        self.rows=default_rows(self.species or None)
        self._refresh()

    def _import_file(self):
        file=filedialog.askopenfilename(parent=self.top,
            title="Import a YAeHMOP parameter file",
            filetypes=[("YAeHMOP parameters","*.txt *.dat *.par *.params"),
                       ("All files","*.*")])
        if not file:return
        try:
            self.rows=file_rows(file)
            self._refresh()
        except Exception as exc:
            messagebox.showerror("Invalid parameter file",str(exc),parent=self.top)

    def _save(self):
        try:
            # Validate before selecting a path, including duplicate shells.
            from .parameter_editor_model import native_table_text
            native_table_text(self.rows)
        except Exception as exc:
            messagebox.showerror("Invalid shell table",str(exc),parent=self.top)
            return
        file=filedialog.asksaveasfilename(parent=self.top,
            title="Save a native YAeHMOP parameter table",
            initialfile="custom_yaehmop_parameters.txt",
            defaultextension=".txt",
            filetypes=[("Text parameter table","*.txt"),("YAeHMOP table","*.dat"),
                       ("All files","*.*")])
        if not file:return
        try:
            path=save_native_table(self.rows,file)
        except Exception as exc:
            messagebox.showerror("Could not save parameters",str(exc),parent=self.top)
            return
        if self.on_save:
            self.on_save(str(path))
        self.top.destroy()


def show_atomic_parameter_editor(parent,species=None,source_file=None,on_save=None):
    return AtomicParameterEditor(parent,species,source_file,on_save)
