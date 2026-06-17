# DELPHI Analysis

Small analysis workspace for the merged SDST, RAWSDST, and RAWFADANA NanoAOD workflow.

## Setup

```bash
micromamba create -y -f environment.yml
source setup.sh
```

## Workflow

The scripts use `SAMPLE_SET = "100kTest"` by default. Edit the few constants at the top of each script when switching samples or files.

```bash
./scripts/merge_nanoaod.py
./scripts/hadd_nanoaod.py
./scripts/plot_variables.py
```

The output layout is:

```text
output/merged/<sample-set>/<sample>/job_<N>.root
output/hadd/<sample-set>/<sample>.root
plots/<sample-set>/
```

`plot_variables.py` takes explicitly labelled input files through `INPUT_SAMPLES`, so plot legends and summaries do not depend on file names.
