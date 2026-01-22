# Jupyter notebook tunnel (cheatsheet)

**Purpose:** View a loom quickly from a remote VM in your local browser, and optionally look inside with Python.

---

## 1) Set up on the VM

```bash
# (optional) pin NumPy if some tools expect < 2.0
python3 -m pip install --user 'numpy<2.0'
python3 -c "import numpy as np; print(np.__version__)"

# install and start Jupyter without opening a browser on the VM
python3 -m pip install --user notebook
jupyter notebook --no-browser --port=8888

# copy the token URL printed by Jupyter
# (you can also list active servers via: jupyter notebook list)
```

## 2) Create an SSH tunnel from your local machine
```bash
# open a new terminal on your local machine
# replace <user>@<host> with your VM login
ssh -N -L 8888:localhost:8888 <user>@<host>

# example (commented; replace with your own details if desired):
# ssh -N -L 8888:localhost:8888 name@virtual.machine.edu
```

## 3) Open the notebook in your local browser
```bash
http://localhost:8888/tree?token=<paste-token-from-remote>
# The token changes each time you start Jupyter—copy it from the VM’s Jupyter output.
```

## 4) Optional: quick Python look at a loom
```python
# Inside a notebook cell (or any Python REPL with loompy installed):
# pip install loompy  # if needed
import loompy

# example: SCENIC input loom in this repo layout
ds = loompy.connect("resources/pyscenic_resources/Joanito_StromaOnly_Filtered_Seurat_Object_All_Cohorts_FullDataSet_SCopeLoomR.loom", mode="r")
list(ds.attrs.items())   # show global attributes
ds.close()

# Minimal REPL example (generic):
>>> import loompy
>>> ds = loompy.connect("cortex.loom")
>>> ds
>>> ds.close()
```


## 5) Notes
- If port 8888 is busy, use 8889 (change both the Jupyter and SSH commands).
- Keep --no-browser on the VM and tunnel only to localhost on your machine.
- Stop Jupyter with Ctrl-C on the VM; stop the SSH tunnel with Ctrl-C locally.
- If your SCENIC-generated loom is not in a “standard” format you expect elsewhere,
  first build the SCENIC input loom with the loom-prep script in this repo:
  scripts/01_loom_prep/01_tas_seurat_filtered_to_scenic_loom.R
