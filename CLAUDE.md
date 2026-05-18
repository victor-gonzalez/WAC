# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

WAC ("Wayne State Analysis Code") — a ROOT/C++17 framework for heavy-ion event
generation, reading, and analysis, focused on two-particle (differential)
correlations and balance functions. Despite the repo name (`CLUSTERMODELWAC`),
the active code is the BEST-sampler / Pythia / UrQMD analysis chain plus cluster
batch tooling; the `Cluster/` directory and the cluster model itself are legacy.

## Building

ROOT (with `$ROOTSYS` set, i.e. `thisroot.sh` sourced) and a Pythia8 install
(`$PYTHIA8`) are prerequisites. Build with:

```
source set-WAC-PATH        # exports WAC_ROOT / WAC_BIN / WAC_LIB, prepends bin/ and lib/ to paths
mkdir build && cd build
cmake ../
make install               # installs libs to ./lib, executables to ./bin
```

Use the site wrapper instead of `set-WAC-PATH` directly when on a known machine —
each one sources ROOT, sets `$PYTHIA8`, then sources `set-WAC-PATH`:
`set-WAC-HOME`, `set-WAC-UCM`, `set-WAC-WSU`, `set-WAC-GSI`. On clusters the
equivalent is `Clusters/<SITE>/entergenwac.sh`.

The top-level `CMakeLists.txt` builds only `Base`, `Best`, `Pythia`, `Urqmd`,
`Utils`, and conditionally `Home`/`GSI` (added when `$MYHOME`/`$ATGSI` is set —
these provide a vendored `EGPythia8`/`TPythia8`). `Pythia` is skipped when
`$ATWSU` is set. The other top-level directories (`Cluster`, `Blastwave`,
`Eccentricity`, `Ampt`, `Music`, `Plots`, `StatStudies`, `CollisionGeometry`,
`Epos`, `LambdaAnalysis`) are **not** wired into the build.

There is no test suite. `build/tests.sh` and `tests/` only exercise the bash
JSON-array helpers in `Clusters/utilities.sh`.

## Running an analysis

Executables live in `bin/` after `make install`. They read a
`configuration.json` from the current working directory (deserialized into a
`BestAnalysisConfiguration` / `PythiaAnalysisConfiguration` via ROOT's
`TBufferJSON`) and write `.root` output under `Output/`. Example:

```
RunBestSimulationTwoParticlesDiff <oscar_input_file> <n_events> <me|nome> <job_index>
```

(`me` enables mixed-event background.) The `Test*AnalysisConfiguration`
executables generate/print a valid `configuration.json`. `statUncertain` /
`statUncertainPythia` / `statUncertainBest` post-process a set of sub-sample
result files into final histograms with statistical uncertainties — see next
section for why this is a separate step.

## Architecture

**Task / EventLoop pipeline.** `Base/Task.hpp` is the unit of work
(`initialize` / `execute` / `finalize`, plus histogram lifecycle hooks). An
`EventLoop` (a `TaskCollection`) holds an ordered list of Tasks and calls
`execute()` on each, once per event, for N events. A run is assembled in a
`Run*Simulation*.cpp` `main()`: build the config, create one *event reader*
Task, create many *analyzer* Tasks, add them all to the `EventLoop`, call
`run(nEvents, nReport)`.

**Single shared Event.** `Base/Event.hpp` is a singleton (`Event::getEvent()`).
The reader Task fills it each iteration; analyzer Tasks read it. Readers:
`Best/BestEventReader` (OSCAR 2013 reduced format), `Pythia/PythiaEventGenerator`,
`Urqmd/UrQMDEvent`, `Ampt/AmptEventReader`. Analyzers (in `Base/`):
`TwoPartDiffCorrelationAnalyzer`, `TwoPartDiffCorrelationAnalyzerME` (mixed
events), `TwoPartCorrelationAnalyzer`, `ParticleAnalyzer`.

**Rapidity vs pseudorapidity is a template parameter.** `ParticleFilter<r>`,
`BestEventReader<r>`, the analyzers, etc. are templated on
`AnalysisConfiguration::RapidityPseudoRapidity` (`kRapidity` / `kPseudorapidity`),
and `TwoPartDiffCorrelationAnalyzer<r, options>` additionally on
`FillPairOptions` (`kNoAdditionalOptions`, `kFillPratt`, `kFillInvariantMass`,
`kFillPrattAndInvariantMass`). `main()` therefore branches into a small tree of
explicit template instantiations — see `Best/RunBestSimulationTwoParticlesDiff.cpp`
for the canonical pattern.

**Filters.** `EventFilter` (`MinBias`, `Centrality` with a percentile range),
`ParticleFilter<r>` (species/charge/pt/rapidity windows), `ParticlePairFilter`.
Note: for balance functions to come out right, particle filters must be ordered
charged-first, with the positive species immediately before its negative
counterpart, neutrals last (documented in `RunBestSimulationTwoParticlesDiff.cpp`).

**Histograms.** `HistogramCollection` / `Histograms` subclasses. `Base` is built
with `-DOPTIMIZEADDBINCONTENT`, which **disables the Sumw2 error structure** for
speed in the pair inner loop — so per-run histogram errors are meaningless;
uncertainties are obtained by running many sub-samples and combining them with
the `statUncertain*` executables.

**ROOT dictionaries.** Each library directory has a `*LinkDef.h` and uses
`ROOT_GENERATE_DICTIONARY`; classes that need to be streamed/JSON-serialized
have `ClassDef`. When adding a streamed class, add it to both `CMakeLists.txt`
(SOURCE/HEADER lists) and the directory's `LinkDef.h`.

## Cluster batch workflow

`Clusters/<SITE>/` (`GSI`, `WSU`, `UCM`) holds SLURM submission scripts.
`batchRun*.sh` is the entry point: it stamps a dated production directory, copies
the `configuration.json`, reads the rapidity / pT-range arrays out of it with the
`extract_json_*` bash functions from `Clusters/utilities.sh`, submits a job array
of `runBest*Correlations.sh` / `runPythiaCorrelations.sh` tasks, then chains
`afterany`-dependent merge jobs (`runMergeResults.sh`, `runMergeSubsamples.sh`)
and finally `batchRunStatsUncertain.sh`. At GSI everything runs inside a CVMFS
Singularity container (`runScriptInSingularity.sh`); the BEST sampler streams
particles to the reader through a named pipe (`runSamplerAndSmash.sh` →
`runBestDataCollectingEngine.sh`). `Best/configuration.json` is the production
config consumed by these scripts. `TOOLS/producing_at_gsi.txt` has the current
manual update/launch checklist.

## Conventions

C++17. Code style is `.clang-format` (Google-based, Linux braces, 2-space indent,
no column limit, includes not sorted). Match the surrounding files: a leading
author/copyright comment block, `WAC_<Name>` include guards, doxygen-ish banner
comments between sections.
