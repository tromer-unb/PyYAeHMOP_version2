"""Isolated compute subprocess for the GUI, preserving numerical engine calls.

The GUI invokes this as a Python module so slow runs never freeze the window
and multiprocessing is not started from the Tk event thread.
"""
import argparse
import json
from pathlib import Path
from .cli import _tolist


def run_gui_calculation(config):
    config=Path(config).expanduser().resolve()
    if config.suffix.lower()==".spin":
        from .SPIN_TEST.cli import run_case
        result,job=run_case(config)
        return job.outdir/"spin_soc_result.json"
    from .friendly import read_easy_job,calculate_easy_job
    from .export import export_csv
    from .plots import export_plots
    job=read_easy_job(config)
    report=calculate_easy_job(job)
    job.output_dir.mkdir(exist_ok=True,parents=True)
    # Keep JSON for GUI visualization even when CSV/PNG-only was requested.
    target=job.output_dir/"summary.json"
    target.write_text(json.dumps(report,default=_tolist,ensure_ascii=False,indent=2)+"\n",
                      encoding="utf-8")
    if "csv" in job.outputs:
        export_csv(report,job.output_dir)
    if "png" in job.outputs:
        export_plots(report,job.output_dir)
    return target


def main(argv=None):
    parser=argparse.ArgumentParser(description="Run a GUI-created .job/.spin input")
    parser.add_argument("config",type=Path)
    args=parser.parse_args(argv)
    try:
        result=run_gui_calculation(args.config)
    except Exception as exc:
        parser.exit(2,f"Calculation failed: {exc}\n")
    print("RESULT_FILE="+str(result),flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
