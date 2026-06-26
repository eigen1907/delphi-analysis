export PROJECT_ROOT=$(realpath $(dirname ${BASH_SOURCE}))
export PYTHONPATH=${PROJECT_ROOT}/python:${PYTHONPATH}
eval "$(micromamba shell hook --shell bash)"
micromamba activate delphi-analysis-py312
