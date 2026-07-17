#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/MACSEtry.py"
ENV_NAME="${MACSE_ENV_NAME:-macse_env}"
MINIFORGE_DIR="${MACSE_MINIFORGE_DIR:-${SCRIPT_DIR}/.miniforge}"

show_help() {
    cat <<'EOF'
Run MACSE multiple sequence alignment for grouped FASTA files.

Required:
  -i, --input DIR      Input directory with FASTA files
  -o, --output DIR     Output directory for MACSE alignments

Optional:
  -e, --env-name NAME  Conda environment name. Default: macse_env
  -h, --help           Show this help message

Environment variables:
  MACSE_ENV_NAME       Default conda environment name
  MACSE_MINIFORGE_DIR  Local Miniforge install directory if conda is absent

Example:
  ./run_macse.sh \
    -i /home/hellstrom/Загрузки/grouped_by_germlines/v \
    -o /home/hellstrom/Документы/MACSE/v

The script finds conda, installs local Miniforge if conda is absent, creates the
conda environment if needed, installs MACSE into it, and runs MACSEtry.py there.
EOF
}

need_cmd() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $1" >&2
        exit 1
    fi
}

download_with_python() {
    local url="$1"
    local output="$2"
    python3 - "$url" "$output" <<'PY'
import sys
import urllib.request

url, output = sys.argv[1], sys.argv[2]
print(f"Downloading {url}", file=sys.stderr)
with urllib.request.urlopen(url) as response, open(output, "wb") as handle:
    handle.write(response.read())
PY
}

find_conda() {
    if command -v conda >/dev/null 2>&1; then
        command -v conda
        return
    fi

    if [[ -x "${MINIFORGE_DIR}/bin/conda" ]]; then
        printf '%s\n' "${MINIFORGE_DIR}/bin/conda"
        return
    fi

    return 1
}

install_miniforge() {
    need_cmd python3

    local os
    local arch
    local installer
    local url
    local installer_path

    os="$(uname -s)"
    arch="$(uname -m)"

    case "${os}:${arch}" in
        Linux:x86_64) installer="Miniforge3-Linux-x86_64.sh" ;;
        Linux:aarch64|Linux:arm64) installer="Miniforge3-Linux-aarch64.sh" ;;
        Darwin:x86_64) installer="Miniforge3-MacOSX-x86_64.sh" ;;
        Darwin:arm64) installer="Miniforge3-MacOSX-arm64.sh" ;;
        *)
            echo "ERROR: unsupported platform for automatic Miniforge install: ${os} ${arch}" >&2
            echo "Install conda manually, then rerun this script." >&2
            exit 1
            ;;
    esac

    url="https://github.com/conda-forge/miniforge/releases/latest/download/${installer}"
    mkdir -p "${SCRIPT_DIR}/.downloads"
    installer_path="${SCRIPT_DIR}/.downloads/${installer}"

    if [[ ! -f "${installer_path}" ]]; then
        download_with_python "${url}" "${installer_path}"
    fi

    echo "Installing Miniforge into ${MINIFORGE_DIR}" >&2
    bash "${installer_path}" -b -p "${MINIFORGE_DIR}"
}

ensure_conda() {
    local conda_bin

    if conda_bin="$(find_conda)"; then
        printf '%s\n' "${conda_bin}"
        return
    fi

    echo "Conda was not found. Installing local Miniforge..." >&2
    install_miniforge
    find_conda
}

conda_env_exists() {
    local conda_bin="$1"
    "${conda_bin}" env list | awk '{print $1}' | grep -qx "${ENV_NAME}"
}

ensure_env() {
    local conda_bin="$1"

    if ! conda_env_exists "${conda_bin}"; then
        echo "Creating conda environment '${ENV_NAME}' with Python and MACSE..."
        "${conda_bin}" create -y -n "${ENV_NAME}" \
            -c conda-forge -c bioconda \
            python=3.12 macse
        return
    fi

    echo "Conda environment '${ENV_NAME}' already exists."
    if ! "${conda_bin}" run -n "${ENV_NAME}" macse -h >/dev/null 2>&1; then
        echo "Installing MACSE into '${ENV_NAME}'..."
        "${conda_bin}" install -y -n "${ENV_NAME}" \
            -c conda-forge -c bioconda \
            macse
    fi
}

input_dir=""
output_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            input_dir="${2:-}"
            shift 2
            ;;
        -o|--output)
            output_dir="${2:-}"
            shift 2
            ;;
        -e|--env-name)
            ENV_NAME="${2:-}"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            show_help >&2
            exit 2
            ;;
    esac
done

if [[ ! -f "${PY_SCRIPT}" ]]; then
    echo "ERROR: MACSEtry.py not found: ${PY_SCRIPT}" >&2
    exit 1
fi

if [[ -z "${input_dir}" ]]; then
    echo "ERROR: input directory is required. Use -i or --input." >&2
    show_help >&2
    exit 2
fi

if [[ -z "${output_dir}" ]]; then
    echo "ERROR: output directory is required. Use -o or --output." >&2
    show_help >&2
    exit 2
fi

conda_bin="$(ensure_conda)"
ensure_env "${conda_bin}"

echo "Running MACSE alignment in conda environment '${ENV_NAME}'..."
echo "Input: ${input_dir}"
echo "Output: ${output_dir}"

"${conda_bin}" run -n "${ENV_NAME}" python "${PY_SCRIPT}" \
    --input "${input_dir}" \
    --output "${output_dir}"
