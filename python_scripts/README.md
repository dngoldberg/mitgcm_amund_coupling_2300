# Coupling Scripts README

This folder contains a collection of shell and SLURM scripts used to build, prepare, run, and manage coupled ocean–ice simulations (likely MITgcm-based) on an HPC system (e.g. ARCHER/Cray environment).

---

## 📁 File Overview

### `pad_with_zeros.sh`
Utility function for formatting numbers.

**What it does:**
- Pads a given number with leading zeros to a specified length.

**Usage:**
```sh
source pad_with_zeros.sh
pad_with_zeros 7 3   # Output: 007
```

---

### `prepare_run.sh`
Prepares clean run directories for both ocean and ice simulations.

**What it does:**
- Creates or resets:
  - `run_ice_<params>`
  - `run_oce_<params>`
- Clears old files and recreates `diags/` directories
- Copies required pickup/restart files into the ocean run directory

**Usage:**
```sh
./prepare_run.sh <arg1> <arg2> <arg3> <arg4> <arg5>
```

---

### `make_archer_oce.sh`
Build script for the ocean model on a Cray/ARCHER system.

**What it does:**
- Loads required modules (Cray environment)
- Creates/cleans `build_oce/`
- Prepares build directory for compilation

---

### `make_archer_ice.sh`
Build script for the ice model.

**What it does:**
- Loads GNU programming environment
- Configures PETSc paths
- Creates/cleans `build_ice/`
- Prepares for compilation

---

### `submit.sh`
Main driver script for launching simulations.

**What it does:**
- Parses input arguments
- Sets defaults for optional parameters
- Prepares run directories (unless restarting)
- Submits the job via SLURM (`run_repeat.slurm`)

---

### `run_repeat.slurm`
SLURM batch script for running the main simulation.

**What it does:**
- Requests HPC resources
- Loads required scientific libraries
- Executes the simulation with provided parameters

---

### `rput_couple_mds.slurm`
Post-processing and data transfer script.

**What it does:**
- Combines model outputs
- Handles synchronization or transfer of results
- Logs timing information

---

## 🔍 Parameter Reference (`submit.sh`)

The `submit.sh` script is the main entry point for launching simulations. It parses command-line arguments and forwards them to the SLURM job (`run_repeat.slurm`) in a fixed positional order.

---

### 🧾 Required Parameters

```bash
./submit.sh -y <year> -c <case> -a <amplitude> -p <project_code>
```

| Flag | Variable     | Meaning (inferred) |
|------|-------------|------------------|
| `-y` | `y_value`   | Year or experiment ID (likely controls forcing or initial conditions) |
| `-c` | `cal_value` | Case / configuration type (e.g. `TC`) |
| `-a` | `PAS_value` | Amplitude or forcing strength parameter |
| `-p` | `p_value`   | Project/account code for SLURM allocation |

---

### ⚙️ Optional Parameters

```bash
[-n <param_name>] [-f <cfric>] [-s <shitrans>] [-d <depth_code>] [-i <iter>] [-r]
```

| Flag | Variable   | Default           | Meaning (inferred) |
|------|-----------|------------------|--------------------|
| `-n` | `parmn`   | `iceParmDefault` | Ice parameter set name |
| `-f` | `cfric`   | `0.01`           | Friction coefficient |
| `-s` | `shitrans`| `0.00014`        | Shelf/ice transfer coefficient |
| `-d` | `depth`   | `G`              | Depth or geometry code |
| `-i` | `iter`    | `0`              | Restart iteration index |
| `-r` | `r_flag`  | `false`          | Restart flag (skip cleanup if set) |

---

### 🔄 Internal Fixed Parameters

The script calls the SLURM job like this:

```bash
run_repeat.slurm \
  $y_value \
  $cal_value \
  $PAS_value \
  coul \
  20 \
  $p_value \
  $cfric \
  $shitrans \
  $depth \
  $parmn \
  $iter
```

| Position | Value   | Meaning (inferred) |
|----------|--------|--------------------|
| 4        | `coul` | Coupling mode (likely ocean–ice coupling enabled) |
| 5        | `20`   | Fixed parameter (possibly coupling interval or timestep block size) |

---

### 🧠 Behaviour Notes

- **Default handling:**  
  Missing optional parameters are automatically assigned default values.

- **Restart mode (`-r`):**
  - Without `-r`: run directories are cleaned and recreated
  - With `-r`: existing run continues without cleanup

- **Formatting:**
  - Numeric values like `.01` are normalized to `0.01`

- **Job naming convention:**
  ```
  c<case><amplitude><param_name>
  ```
  Example:
  ```
  cTC20iceParm4
  ```

---

### ✅ Example Usage

```bash
./submit.sh \
  -y 2009 \
  -c TC \
  -a 20 \
  -p n02-GRISLAKES \
  -n iceParm4 \
  -f 0.014 \
  -s 0.00014 \
  -d G \
  -i 25 \
  -r
```

---

## 🔄 Typical Workflow

1. Build models:
```sh
./make_archer_oce.sh
./make_archer_ice.sh
```

2. Prepare run directories:
```sh
./prepare_run.sh <params>
```

3. Submit simulation:
```sh
./submit.sh [options]
```

4. Run simulation via SLURM:
```sh
sbatch run_repeat.slurm
```

5. Post-process results:
```sh
sbatch rput_couple_mds.slurm
```

---

## ⚙️ Environment Requirements

- HPC system (Cray/ARCHER-like)
- SLURM scheduler
- Required modules:
  - `PrgEnv-cray` or `PrgEnv-gnu`
  - `cray-hdf5-parallel`
  - `cray-netcdf-hdf5parallel`
- PETSc (for ice model)

---

## ⚠️ Notes

- Scripts rely on a fixed directory structure (`run_*`, `build_*`)
- Environment variables (e.g. `ROOTDIR`) must be set correctly
- Input/restart files must exist before running
- Parameter meanings are inferred and should be validated against model configs

---

## 🛠️ Tips

- Test small runs before scaling up
- Keep logs for debugging
- Verify module versions on your HPC system
- Inspect `run_repeat.slurm` for deeper parameter usage
