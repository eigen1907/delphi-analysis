export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${PROJECT_ROOT}/python${PYTHONPATH:+:${PYTHONPATH}}"
uv sync --locked --project "${PROJECT_ROOT}" || return 1
source "${PROJECT_ROOT}/.venv/bin/activate"
