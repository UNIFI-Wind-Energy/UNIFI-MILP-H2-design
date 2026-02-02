# MILP Microgrid Optimization (Wind–PV–BESS–Hydrogen) with Gurobi

This repository contains a **mixed-integer linear programming (MILP)** model to optimize the **operation** and optionally the **design sizing** of a hybrid energy system including:

- Wind power (WT)
- Photovoltaics (PV)
- Battery energy storage system (BESS)
- Hydrogen chain (electrolyzer + H2 tank + fuel cell)
- Grid exchange (import/export)

The objective is to **maximize the Net Present Value (NPV)** over a multi-year project lifetime using hourly time series for load and renewable generation.

---

## Citation notice

If you use the **MILP optimization framework**, please cite:

F. Superchi and A. Bianchini, *Renewable Energy*, vol. 256, p. 124470, Jan. 2026.  
DOI: 10.1016/j.renene.2025.124470

If you use the **market analysis results and component cost assumptions** implemented in this code, please cite:

F. Superchi, A. Bianchini, A. Moustakis, and G. Pechlivanoglou, *Renewable Energy*, vol. 245, p. 122813, Jun. 2025.  
DOI: 10.1016/j.renene.2025.122813

---

## Disclaimer

- This code is provided for research and educational use.
- The authors assume no responsibility for results obtained from using or modifying the code.
- Users are responsible for validating assumptions, inputs, and outputs for their specific application.

---

## Requirements

### Python packages

- `numpy`
- `pandas`
- `gurobipy`

Install the Python dependencies (example):

```bash
pip install numpy pandas
```

### Gurobi

This project requires **Gurobi Optimizer** (commercial solver) and a valid license.

- Gurobi is **not distributed** with this repository.
- Users must install Gurobi separately and set up a valid license.

Quick check that `gurobipy` is available:

```bash
python -c "import gurobipy as gp; print(gp.gurobi.version())"
```

---

## Inputs

### 1) Hourly time series (mandatory)

The script requires **three hourly time series** covering the same period:

- `Load.csv`  (electrical demand)
- `PV.csv`    (PV electrical power)
- `WT.csv`    (wind electrical power)

Expected characteristics:

- Resolution: **1 hour**
- Same timestamps for all series
- Same unit for all power series (for example kW)

How timestamps are handled in the current script:

- The timestamps are created internally using `pd.date_range(...)`.
- The CSV files are expected to contain **only the values** (one column each), aligned with the generated hourly index.

If the CSV files include timestamps, adapt the loading block accordingly (recommended for robustness).

### 2) Component costs (mandatory)

The script reads component CAPEX values from:

- `prices_excel.xlsx`

Sheets expected (as currently coded):

- `Li-BESS`
- `ALK EL`
- `PEM FC`
- `H2 Tank`
- `PV`
- `Onshore WT`

Columns expected:

- `avg` (used in the code)

The variable `year` selects which row of the table is used:

- `year = 2020`  -> first row (`k = 0`)
- `year = 2030`  -> second row (`k = 1`)
- `year = 2050`  -> third row (`k = 2`)

---

## Model modes

Two switches control the behavior:

- `MILP_design = True`  
  Design sizing variables are optimized (BESS capacity, electrolyzer power, fuel cell power, tank capacity, PV scale).

- `MILP_design = False`  
  Fixed design is used (user-defined values in the script).

The script currently uses:

```python
MILP_design = True
kWh_factor = 1  # hourly simulation
```

---

## Key assumptions (important)

- Hourly timestep (`kWh_factor = 1`).
- No simultaneous charge and discharge for the battery (binary constraint).
- No simultaneous hydrogen production and consumption (binary constraint).
- State-of-charge boundary conditions are enforced:
  - Battery starts at 50% SOC
  - H2 tank starts and ends at 40% SOC
- Grid import and export are limited by `P_grid_max`.
- Export remuneration uses band and over-band power and includes remuneration from BESS export.

---

## Running the code

1) Put the following files in the working directory:

- `Load.csv`
- `PV.csv`
- `WT.csv`
- `prices_excel.xlsx`

2) Run:

```bash
python main.py
```

---

## Outputs

After a successful optimization, the script writes:

- `hourly.csv`  
  Hourly dispatch and state trajectories (production, SOC, grid, curtailment, etc.)

- `sizes.csv`  
  Optimal component sizes (only if `MILP_design = True`)

- `NPV.csv`  
  Discounted cash flow terms per year (as implemented in the script)

- `results.csv`  
  Aggregated annual energies and revenues

All outputs are saved as CSV using `;` as separator.

---

## Troubleshooting

### `ImportError: No module named gurobipy`

Gurobi (and its Python API) is not installed or not visible in the current Python environment.

- Install Gurobi and enable `gurobipy` in the environment used to run the script.
- Verify with:

```bash
python -c "import gurobipy as gp; print(gp.gurobi.version())"
```

### Gurobi license errors

A valid Gurobi license must be installed and configured on the machine.

### Infeasible model

If the solver returns infeasible or unbounded:

- Verify that input series contain numeric values and no missing rows.
- Check unit consistency across `load`, `PV`, and `WT`.
- Check that sizing bounds and `P_grid_max` are compatible with the load magnitude.
- Temporarily relax end-of-horizon SOC constraints to diagnose infeasibility.

---

## Suggested next improvements (optional)

- Replace the value-only CSV format with timestamped CSV parsing and automatic alignment.
- Move parameters into a YAML or JSON configuration file.
- If solver-agnostic execution is required, add a modeling layer that can target open-source solvers (for example HiGHS or CBC).
