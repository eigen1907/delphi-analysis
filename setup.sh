export PROJECT_PREFIX=$(realpath $(dirname ${BASH_SOURCE}))
export PYTHONPATH=${PROJECT_PREFIX}/python:${PYTHONPATH}
eval "$(micromamba shell hook --shell bash)"
micromamba activate delphi-analysis
