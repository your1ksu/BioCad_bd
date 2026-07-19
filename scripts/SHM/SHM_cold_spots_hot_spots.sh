#!/usr/bin/env bash

# Завершать работу при любой ошибке
set -e

# --- НАСТРОЙКИ CONDA ---
ENV_NAME="shm_env"
PYTHON_VERSION="3.11"
MINICONDA_URL="https://anaconda.com"
CONDA_DIR="$HOME/miniconda3"

# Путь к самому Python-скрипту по умолчанию
SCRIPT_PATH="$HOME/Документы/SHM/SHM_cold_spots_hot_spots.py"

# --- ФУНКЦИЯ СПРАВКИ ---
show_help() {
    echo "Использование: $0 -i <путь_к_fasta> -o <путь_к_json> [-s <путь_к_python_скрипту>]"
    echo ""
    echo "Параметры:"
    echo "  -i, --input    Путь к входному FASTA-файлу или папке с файлами"
    echo "  -o, --output   Путь к выходному JSON-файлу отчета"
    echo "  -s, --script   Путь к Python-скрипту анализа"
    echo "  -h, --help     Показать это справочное сообщение"
    exit 0
}

# --- ПАРСИНГ КЛЮЧЕЙ КОМАНДНОЙ СТРОКИ ---
INPUT_PATH=""
OUTPUT_PATH=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT_PATH="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_PATH="$2"
            shift 2
            ;;
        -s|--script)
            SCRIPT_PATH="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        *)
            echo "Ошибка: Неизвестный параметр $1"
            show_help
            exit 1
            ;;
    esac
done

# Проверяем, переданы ли обязательные ключи
if [ -z "$INPUT_PATH" ] || [ -z "$OUTPUT_PATH" ]; then
    echo "Ошибка: Ключи --input и --output являются обязательными."
    show_help
    exit 1
fi

echo "=== Проверка окружения для скрипта SHM ==="

# 1. Проверяем, установлена ли Conda
if ! command -v conda &> /dev/null; then
    if [ -f "$CONDA_DIR/bin/conda" ]; then
        echo "Conda найдена в директории: $CONDA_DIR"
        source "$CONDA_DIR/etc/profile.d/conda.sh"
    else
        echo "Conda не найдена. Начинается загрузка и установка Miniconda..."
        TMP_SH=$(mktemp)
        curl -sL "$MINICONDA_URL" -o "$TMP_SH"
        bash "$TMP_SH" -b -p "$CONDA_DIR"
        rm "$TMP_SH"
        source "$CONDA_DIR/etc/profile.d/conda.sh"
        conda init bash
        echo "Miniconda успешно установлена!"
    fi
else
    echo "Conda уже установлена в системе."
    eval "$(conda shell.bash hook)"
fi

# 2. Проверяем и создаем виртуальную среду
if conda info --envs | grep -q "$ENV_NAME"; then
    echo "Среда Conda '$ENV_NAME' уже существует."
else
    echo "Создание новой среды Conda '$ENV_NAME' с Python $PYTHON_VERSION..."
    conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y
fi

# 3. Активируем среду
echo "Активация среды '$ENV_NAME'..."
conda activate "$ENV_NAME"

# 4. Проверяем наличие Python-скрипта и запускаем анализ с переданными ключами
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Ошибка: Файл скрипта не найден по пути $SCRIPT_PATH"
    exit 1
fi

echo "Запуск анализа SHM..."
python3 "$SCRIPT_PATH" --input "$INPUT_PATH" --output "$OUTPUT_PATH"

echo "=== Процесс успешно завершен! ==="
