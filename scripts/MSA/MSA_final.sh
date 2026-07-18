#!/usr/bin/env bash
set -euo pipefail

# Bootstrap and run ../MSA/MSA_final.py in a conda environment with MAFFT.
#
# Usage:
#   bash MSA.sh -i /path/to/input_fastas -o /path/to/aligned_output
#
# The wrapped Python script uses only the Python standard library. Its only
# external runtime dependency is the MAFFT executable.

ENV_NAME="${MSA_ENV_NAME:-msa_final_env}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MSA_PY="${MSA_PY:-${SCRIPT_DIR}/../MSA/MSA_final.py}"
MINIFORGE_DIR="${MSA_MINIFORGE_DIR:-${SCRIPT_DIR}/.miniforge}"

if [[ ! -f "${MSA_PY}" ]]; then
    echo "ERROR: MSA_final.py not found: ${MSA_PY}" >&2
    echo "Set MSA_PY=/path/to/MSA_final.py or place this script next to ../MSA/MSA_final.py." >&2
    exit 1
fi

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
    local installer_path="${SCRIPT_DIR}/.downloads/${installer}"

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

ensure_env() {
    local conda_bin="$1"

    if "${conda_bin}" env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
        return
    fi

    echo "Creating conda environment '${ENV_NAME}' with Python and MAFFT..."
    "${conda_bin}" create -y -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        python=3.12 mafft
}

main() {
    local conda_bin
    conda_bin="$(ensure_conda)"
    ensure_env "${conda_bin}"

    echo "Running MSA_final.py in environment '${ENV_NAME}'..."
    "${conda_bin}" run -n "${ENV_NAME}" python "${MSA_PY}" "$@"
}

main "$@"
