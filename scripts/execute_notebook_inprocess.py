"""Execute notebook code cells in order without a networked Jupyter kernel.

Useful on runners that do not allow local kernel sockets. Uses real IPython
execution and captures actual rich outputs; never synthesizes analysis output.
This does not verify the Jupyter server/kernel transport.
"""
import argparse
from datetime import datetime, timezone
from pathlib import Path

import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("notebook", nargs="?", default="notebooks/01_pbmc3k_walkthrough.ipynb")
    args = p.parse_args()
    path = Path(args.notebook).resolve()
    nb = nbformat.read(path, as_version=4)
    shell = InteractiveShell.instance()
    count = 0
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        count += 1
        print(f"Execute code cell {count} (notebook cell {i + 1})", flush=True)
        with capture_output(stdout=True, stderr=True, display=True) as captured:
            result = shell.run_cell(cell.source, store_history=True)
        if result.error_before_exec or result.error_in_exec:
            raise RuntimeError(f"Cell {count} failed: {captured.stdout}\n{captured.stderr}") from (result.error_before_exec or result.error_in_exec)
        outputs = []
        if captured.stdout:
            outputs.append(nbformat.v4.new_output("stream", name="stdout", text=captured.stdout))
        if captured.stderr:
            outputs.append(nbformat.v4.new_output("stream", name="stderr", text=captured.stderr))
        for item in captured.outputs:
            outputs.append(nbformat.v4.new_output("display_data", data=item.data, metadata=item.metadata))
        cell.outputs = outputs
        cell.execution_count = count
    nb.metadata["validation"] = {"mode": "real code execution via in-process IPython; Jupyter socket transport not tested",
                                  "completed_at_utc": datetime.now(timezone.utc).isoformat(), "executed_code_cells": count}
    nbformat.validate(nb)
    nbformat.write(nb, path)
    print(f"Successfully executed {count} code cells and validated notebook format.")


if __name__ == "__main__":
    main()
