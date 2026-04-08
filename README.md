# GeneTransfer

Simulation code accompanying the manuscript:

> **A minimal model of panimmunity maintenance by horizontal gene transfer in the ecological dynamics of bacteria and phages**
> Wenping Cui, Jemma M. Fendley, Sriram Srikant, Boris I. Shraiman
> *PNAS* 122 (31), e2417628122 (2025)
> DOI: [10.1073/pnas.2417628122](https://doi.org/10.1073/pnas.2417628122)

## Overview

Microbial pan-genomes harbour remarkable diversity of defence and counter-defence genes, shaped by the perpetual arms race between bacteria and bacteriophages.  While no single genome can carry defences against all phages, the **panimmunity hypothesis** posits that bacteria can acquire the necessary defence genes from the surrounding environment through horizontal gene transfer (HGT).

This repository implements a **stochastic population-genetics simulation** of phage–host co-evolution that tests this hypothesis.  The model shows that HGT, above a critical rate threshold, enables the stable maintenance of a diverse set of defence/counter-defence genes in the community pan-genome — even under continuous turnover of individual bacterial and phage strains driven by boom/bust cycles.  The mechanism is formally analogous to MacArthur–Wilson island biogeography: genetic elements "migrate" between bacterial genomes in the same way species migrate between islands.

## Repository structure

```
GeneTransfer/
├── main.jl                      # CLI entry point — parse args, run simulation, save output
├── config.json                  # Default simulation parameters (overridden by CLI flags)
├── src/
│   └── simulation_lightning.jl  # Core simulation library (module simulation_light)
├── test/
│   ├── runtests.jl              # Unit tests (219 tests across all exported functions)
│   └── compare_versions.jl     # Output equivalence check between refactored and original code
├── run_parallel.sh              # Batch parameter sweep using GNU parallel
├── SinglePair_simulation.ipynb  # Single host–phage pair: τ-leap stochastic + LV deterministic
├── ManyPairs_Simulation.ipynb   # Many genotype pairs: full population dynamics & statistics
├── Paper_plots.ipynb            # Publication figures from pre-computed data in Data4Plot/
└── Data4Plot/                   # Pre-computed simulation outputs used for paper figures
```

## Model description

Each organism carries **d** genes drawn from a pool of **L** available genes.  Genotypes are represented as sorted d-element subsets of {1, …, L}.  The simulation evolves host and phage populations forward in discrete generations, applying (in order each generation):

1. **Recombination / HGT** — one or more of the following mechanisms (controlled by rate parameters):
   - Host random mutation (`rR_random`)
   - Host pairwise recombination (`rR_hosthost`)
   - Phage random mutation (`rH_random`)
   - Phage–phage recombination (`rH_phagephage`)
   - Phage–host HGT — 2-body (`rH_phagehost`)
   - Phage–host HGT — 3-body coinfection (`rH_3body`)

2. **Fitness & selection** — Lotka-Volterra-like frequency-dependent fitness; Poisson-sampled offspring under hard, soft, or unconstrained population-size regimes.

**Summary statistics** computed from time series:
- θ diversity at genotype and gene level (host and phage)
- Genotype persistence times and gap times
- Establishment sizes

## Requirements

- [Julia](https://julialang.org/) ≥ 1.9
- Julia packages: `ArgParse`, `BenchmarkTools`, `Combinatorics`, `DataStructures`, `Distributions`, `JSON`, `SpecialFunctions`, `StatsBase`
- For notebooks: `IJulia`, `PyPlot`, `DifferentialEquations`, `ImageFiltering`, `LaTeXStrings`, `OffsetArrays`
- For batch runs: [GNU parallel](https://www.gnu.org/software/parallel/)

## Usage

### Single run

```bash
julia main.jl \
  --T 200000 --L 40 --J 0.005 --N 1E6 \
  --rR_hosthost 1E-3 --rH_phagehost 1E-3 \
  --population_constraint hard \
  --save_dir results/
```

All parameters default to the values in `config.json` and can be overridden on the command line.

### Key parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--L` | Total number of genes in the pool | 30 |
| `--d` | Genes per genotype | 2 |
| `--N` | Host population size | 1 × 10⁶ |
| `--J` | Interaction strength (host–phage coupling) | 5.0 |
| `--T` | Simulation generations | 5 × 10⁵ |
| `--rR_hosthost` | Host–host recombination rate | 0 |
| `--rH_phagehost` | Phage–host HGT rate (2-body) | 0 |
| `--rH_3body` | Phage–host HGT rate (3-body coinfection) | 0 |
| `--population_constraint` | `hard`, `soft`, or `no` | `soft` |
| `--save_dir` | Output directory | `results_simulation_test/` |

### Batch parameter sweep

```bash
bash run_parallel.sh
```

Runs 28 parallel jobs sweeping recombination rates over [10⁻⁶, 10⁻¹·⁵], N ∈ {10⁶, 10⁷}, and L ∈ {30, 40, 50}.

### Merge output files

```bash
awk '!/^(L,)/' results/*.txt > combined_output.txt
```

### Run tests

```bash
julia test/runtests.jl
```

## Contact

For questions about the code, please contact [wenpingcui@gmail.com](mailto:wenpingcui@gmail.com).
