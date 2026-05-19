# Hilbert Basis Algorithm

Code for *Scalable Enumeration of Pareto-optimal Polymers for Computing Equilibrium Concentrations*.

The normal workflow is:

1. write or choose a monomer file;
2. run `src/hilbert_pipeline.py` to enumerate polymer vectors;
3. optionally run `src/coffee_pipeline.py` to compute equilibrium concentrations and export them as CSV.
---

## Repository layout

```text
example-tbns/             Example monomer files and the TBN generator
src/hilbert_pipeline.py   Main Hilbert-basis / covering-design pipeline
src/coffee_pipeline.py    Build COFFEE inputs, run COFFEE, write sorted CSV output
src/coffee_parser.py      Shared COFFEE helpers
src/benchmarks/           Paper benchmark and leakage-analysis scripts
coffee/                   COFFEE source checkout
results/                  Generated outputs; safe to delete/regenerate
```

---

## 1. Set up Python

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install requests beautifulsoup4 lxml numpy matplotlib
```

---

## 2. Install Normaliz

Normaliz is required for the Hilbert-basis computation. It is **not redistributed** in this repository.

Official Normaliz downloads are at:

```text
https://github.com/Normaliz/Normaliz/releases
```

The Normaliz documentation says each release includes executables for Linux 64-bit, macOS, and Windows 64-bit. Version 3.11.1 is the version used for this project.

### Option A: use a prebuilt Normaliz binary

1. Open the Normaliz releases page:

   ```text
   https://github.com/Normaliz/Normaliz/releases
   ```

2. Download the archive for your operating system, for example one of:

   ```text
   normaliz-3.11.1-MacOS.zip
   normaliz-3.11.1-Linux64.zip
   normaliz-3.11.1-Win64.zip
   ```

   The exact asset name may vary slightly by release.

3. Unzip it somewhere outside this repo or under the gitignored path `src/Normaliz/`.

4. Tell the pipeline where the executable is:

   ```bash
   export NORMALIZ_EXE=/absolute/path/to/normaliz
   ```

5. Check it works:

   ```bash
   "$NORMALIZ_EXE" --version
   ```

### Option B: build Normaliz locally in this repo

```bash
git clone https://github.com/Normaliz/Normaliz.git src/Normaliz
cd src/Normaliz
./bootstrap.sh
./configure
make -j
cd ../..
```

Then check:

```bash
src/Normaliz/source/normaliz --version
```

The pipeline looks for Normaliz in this order:

1. `$NORMALIZ_EXE`
2. `src/Normaliz/source/normaliz`
3. `normaliz` on your shell `PATH`

---

## 3. Build COFFEE, optional

COFFEE is only needed if you want equilibrium concentrations after computing polymer vectors. It is included as a git submodule. If `coffee/` is empty, initialize it first:

```bash
git submodule update --init --recursive
```

COFFEE requires Rust/Cargo. If `cargo --version` fails, install Rust first:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

Then build COFFEE:

```bash
cd coffee/crates/coffee-cli
cargo build --release
cd ../../..
```

The executable should be:

```text
coffee/crates/coffee-cli/target/release/coffee-cli
```

Check it works:

```bash
coffee/crates/coffee-cli/target/release/coffee-cli --help
```

---

## Monomer input format

A monomer file has one monomer per line. Domains are whitespace-separated; `*` marks complementary domains.

```text
# comments are allowed
x1_1 x1_2
a1* a2* b1* b2*
a1 a2 b1 b2 c1
```

If a line contains a label followed by `:`, the label is ignored:

```text
input_1: x1_1 x1_2
```

Example monomer files are in `example-tbns/`.

---

## Run the Hilbert pipeline

Basic example:

```bash
cd src
python hilbert_pipeline.py \
  --monomer-file ../example-tbns/monomers_cascade_n7.txt \
  --t 5 \
  --save
```

This probes candidate block sizes `k`, chooses the best projected runtime, runs the covering-design enumeration, and saves polymer vectors under `results/hilbert_output/`. The probe step downloads covering designs from the La Jolla Covering Repository when available, so internet access is recommended.

### Flags

| Flag | Required? | Default | Description |
|---|---:|---|---|
| `--monomer-file PATH` | no | `example-tbns/monomers_cascade_n8.txt` | Input monomer file. Use repo-relative or absolute paths. |
| `--t INT` | no | `5` | Support bound. The covering design guarantees coverage of all polymer candidates using at most `t` monomer types. |
| `--k INT` | no | auto | Fixed covering block size. If omitted, the pipeline probes `k` values and chooses the best one. |
| `--mode monomer` | no | `monomer` | Cover subsets of monomer types. This is the main mode used in the paper. |
| `--mode domain` | no | `monomer` | Cover subsets of domain types instead. Experimental / appendix mode. |
| `--include-base` | no | off | Also include the full-system case `k = n`. Ignored if `--k` is given. |
| `--fallback-dp` | no | off | If LJCR has no direct `(n,k,t)` design, use the bundled GPK dynamic-programming construction data to build one. Requires internet access as well. |
| `--probe` | no | off | Probe runtimes only; do not run the full enumeration. Ignored if `--k` is given. |
| `--save` | no | off | Save polymer-vector output to disk. |
| `--save-dir PATH` | no | `results/hilbert_output/` | Directory for saved polymer-vector files. Used only with `--save`. |

### Usage examples

Probe and run automatically chosen `k`:

```bash
cd src
python hilbert_pipeline.py --monomer-file ../example-tbns/monomers_cascade_n7.txt --t 5 --save
```

Probe only, useful before launching a long run:

```bash
cd src
python hilbert_pipeline.py --monomer-file ../example-tbns/monomers_cascade_n7.txt --t 5 --probe
```

Run a fixed `k` without probing:

```bash
cd src
python hilbert_pipeline.py --monomer-file ../example-tbns/monomers_cascade_n7.txt --t 5 --k 25 --save
```

Save to a specific directory:

```bash
cd src
python hilbert_pipeline.py \
  --monomer-file ../example-tbns/monomers_cascade_n7.txt \
  --t 5 \
  --k 25 \
  --save \
  --save-dir ../results/my_run
```

Include the full Hilbert-basis baseline in the probe sweep:

```bash
cd src
python hilbert_pipeline.py \
  --monomer-file ../example-tbns/monomers_cascade_n7.txt \
  --t 5 \
  --include-base \
  --save
```

Use local DP fallback if an online covering design is missing:

```bash
cd src
python hilbert_pipeline.py \
  --monomer-file ../example-tbns/monomers_cascade_n7.txt \
  --t 5 \
  --k 25 \
  --fallback-dp \
  --save
```

### Covering-design DP fallback

For larger block sizes, LJCR may not store the exact covering design requested by the pipeline. Passing `--fallback-dp` lets the code build the missing covering using the GPK dynamic-programming construction data bundled with this repository:

```text
data/covering_design/gpk_dp.sqlite
```

Use it like this:

```bash
python src/hilbert_pipeline.py \
  --monomer-file example-tbns/monomers_cascade_n7.txt \
  --t 5 \
  --k 30 \
  --fallback-dp \
  --save
```

`--fallback-dp` requires internet access as well because the construction may need to retrieve smaller covering designs from LJCR.

#### Rebuild the bundled construction data, advanced

The bundled database covers `v ≤ 150`, `k ≤ 80`, `t ≤ 8`. To rebuild or extend it, run from the repository root:

```bash
python src/covering_design/rebuild_dp_db.py \
  --bounds-only \
  --out-dir data/covering_design_rebuild \
  --V 150 \
  --K 80 \
  --T 8 \
  --replace
```

This creates a new `gpk_dp.sqlite`; replace `data/covering_design/gpk_dp.sqlite` only after validating the rebuilt database.

Output files contain one polymer vector per line:

```text
0 0 1 0 0 ...
1 0 0 2 0 ...
```

Entry `i` is the multiplicity of monomer `i` in that polymer.

---

## Run COFFEE on a polymer file

After saving polymer vectors, run:

```bash
python src/coffee_pipeline.py \
  --monomers example-tbns/monomers_cascade_n7.txt \
  --polymers results/hilbert_output/<your_polymer_file>.txt \
  --out-dir results/common/coffee/cascade_n7 \
  --coffee-cli coffee/crates/coffee-cli/target/release/coffee-cli \
  --label cascade_n7
```

This generates COFFEE inputs, runs `coffee-cli`, keeps the raw COFFEE output, and also writes a sorted CSV.

### Flags

| Flag | Required? | Default | Description |
|---|---:|---|---|
| `--monomers PATH` | yes | none | The same monomer file used to generate the polymer vectors. |
| `--polymers PATH` | yes | none | Polymer-vector file produced by `hilbert_pipeline.py` or another benchmark script. |
| `--out-dir PATH` | yes | none | Directory where COFFEE input/output files will be written. |
| `--coffee-cli PATH` | yes | none | Path to the compiled COFFEE executable. Usually `coffee/crates/coffee-cli/target/release/coffee-cli`. |
| `--label TEXT` | no | `coffee` | Label used in log messages only. |

### COFFEE outputs

The output directory contains:

```text
input.ocx                  COFFEE polymer/energy input
input.con                  starting monomer concentrations
domain_energies.txt        domain energy table used for the run
coffee_output.txt          raw COFFEE output
coffee_output_sorted.csv   concentrations paired with polymer vectors, sorted high to low
```

`coffee_output_sorted.csv` has this format:

```csv
concentration_M,polymer_vector
2.630000e-09,0 1 1 1 1 1 0 0
```

---

## Reproduce the paper experiments

Use the consolidated benchmark runner from the repository root:

```bash
./run_all_benchmarks.sh
```

By default, verbose benchmark logs are off to avoid very large log files. To keep logs, run:

```bash
./run_all_benchmarks.sh --logs
```

or call the Python runner directly:

```bash
python src/benchmarks/run_benchmarks.py --experiments all
python src/benchmarks/run_benchmarks.py --experiments all --logs   # keep verbose logs
```

You can also run selected experiments:

```bash
python src/benchmarks/run_benchmarks.py --experiments 1              # runtime vs k
python src/benchmarks/run_benchmarks.py --experiments 2              # runtime vs t + Full HB
python src/benchmarks/run_benchmarks.py --experiments 3              # equilibrium recovery
python src/benchmarks/run_benchmarks.py --experiments 4              # leakage analysis
python src/benchmarks/run_benchmarks.py --experiments 1 3            # multiple selected experiments
```

The four experiments are:

1. **Runtime vs k** for `linear_cascade_n7` and `damien_n10`. Cascade-7 uses the DP covering fallback and, after `k=25`, tests `k` in increments of 5.
2. **Runtime vs t** for linear cascades `m=5..9`, binary trees `d=3,4`, and DNA cascades `m=4..7`. For each `(system,t)`, the script probes `k`, records the probe time, runs the best covering enumeration only if the projected time is below 30 minutes, and also runs the Full-HB baseline.
3. **Equilibrium recovery** for cascade-7 with the first input removed: full P* versus `t=3` and `t=5` at `k=25`, using COFFEE at `1 µM` initial concentration and `-20` binding energy.
4. **Leakage analysis**: first compare removed inputs `K=1..7`, then compare leakage for different `t` values against full P*. Only polymers above `1 nM` are counted in leakage summaries.

Outputs are written to:

```text
results/benchmarks/01_runtime_vs_k/
  runtime_vs_k.csv, runtime_vs_k.jsonl, figures/PNG files
results/benchmarks/02_runtime_by_t/
  probe_details.csv, runtime CSVs, runtime_by_t_*.png
results/benchmarks/03_equilibrium_recovery/
  figures/        equilibrium relative-error plots
  csv/            compact output index
  hilbert_basis/  full and reduced polymer-vector files
  coffee/         COFFEE inputs, raw outputs, and sorted polymer CSVs
results/benchmarks/04_leakage/
  removed_inputs/ removed-input sweep summary, plot, per-K CSVs
  vary_t/         t-sweep summary, plot, per-polymer CSVs
  csv/            tidy aggregate CSV exports
```

Shared caches used across experiments live at:

```text
results/common/hilbert_basis/         reusable Hilbert-basis polymer-vector files
results/common/coffee/                reusable COFFEE inputs, raw outputs, and sorted CSVs
results/benchmarks/logs/              verbose logs, only when run with --logs
```

---

## Runtime safeguards

- Each Normaliz call is capped at 3 hours in `hilbert_pipeline.py`.
- Phase 2 benchmark full enumerations are skipped when the probe-projected time is over 30 minutes.
- `normaliz_watchdog.sh` is available as an external backup watchdog.
