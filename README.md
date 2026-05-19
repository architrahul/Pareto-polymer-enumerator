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
| `--fallback-dp` | no | off | If the La Jolla Covering Repository has no `(n,k,t)` design, construct one locally using the GPK dynamic-programming construction. |
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
  --out-dir results/coffee/cascade_n7 \
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

Benchmark scripts live in `src/benchmarks/`. They are separate from the general-use pipeline.

From `src/benchmarks/`:

```bash
# Figures 2, 3, 4 benchmark phases
python experiments.py --phases 1
python experiments.py --phases 2
python experiments.py --phases 3
python make_plots.py

# Figures 5, 6 leakage analysis for cascade n=7 incomplete
python leakage_compute_all.py --n 7
python leakage_coffee.py --n 7
python plot_leakage_figures.py --n 7
```

Key result locations:

```text
results/experiments/                  benchmark JSON and Hilbert bases
results/figures/                      generated paper figures
results/leakage/hilbert_basis/        leakage Hilbert-basis files
results/leakage/coffee/               leakage COFFEE inputs/outputs/CSVs
results/leakage/analysis/             leakage summaries and plots
results/csv/                          tidy CSV exports
```

---

## Runtime safeguards

- Each Normaliz call is capped at 3 hours in `hilbert_pipeline.py`.
- Phase 2 benchmark full enumerations are skipped when the probe-projected time is over 30 minutes.
- `normaliz_watchdog.sh` is available as an external backup watchdog.
