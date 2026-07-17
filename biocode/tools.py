"""Обнаружение и безопасный запуск внешних бинарей (MAFFT/IQ-TREE/MrBayes/IgBLAST).

Правила (см. PHYLO_ARCHITECTURE.md §2):
  • поиск: явный путь из конфига → $PATH → известные conda-env;
  • не найдено → :class:`ToolNotFoundError` с командой установки (не тихий фолбэк);
  • запуск через список аргументов, БЕЗ ``shell=True`` (кросс-платформенность,
    безопасность), с таймаутом и захватом stdout/stderr;
  • версия инструмента логируется в манифест (воспроизводимость).
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from .errors import ToolNotFoundError, ToolRunError
from .logging_ import get_logger

log = get_logger("tools")

# известные conda-env, где команда собрала инструменты
_CONDA_ENV_BINS = [
    Path("/opt/miniconda3/envs/bv2026_msa/bin"),
    Path("/opt/miniconda3/envs/bv2026/bin"),
]

_INSTALL_HINTS = {
    "mafft": "conda install -c bioconda mafft",
    "iqtree": "conda install -c bioconda iqtree",
    "mb": "conda install -c bioconda mrbayes",
    "igblastn": "conda install -c bioconda igblast",
}

# альтернативные имена бинаря
_ALIASES = {
    "iqtree": ["iqtree", "iqtree2", "iqtree3"],
    "mb": ["mb", "mrbayes"],
}


def find_tool(name: str, explicit: str = "") -> Path:
    """Найти бинарь по имени. ``explicit`` — путь из конфига (высший приоритет)."""
    if explicit:
        p = Path(explicit)
        if p.is_file() and os.access(p, os.X_OK):
            return p
        raise ToolNotFoundError(name, f"указанный путь не исполняем: {explicit}")
    for cand in _ALIASES.get(name, [name]):
        which = shutil.which(cand)
        if which:
            return Path(which)
        for envbin in _CONDA_ENV_BINS:
            p = envbin / cand
            if p.is_file() and os.access(p, os.X_OK):
                return p
    raise ToolNotFoundError(name, _INSTALL_HINTS.get(name, ""))


def tool_version(binpath: Path, name: str) -> str:
    """Строка версии инструмента (best-effort, не падаем если не распарсилось)."""
    for flag in ("--version", "-version", "-v"):
        try:
            r = subprocess.run([str(binpath), flag], capture_output=True,
                               text=True, timeout=30)
            out = (r.stdout or "") + (r.stderr or "")
            for line in out.splitlines():
                if any(k in line.lower() for k in ("version", name.lower(), "v1", "v2")):
                    return line.strip()
            if out.strip():
                return out.strip().splitlines()[0]
        except Exception:
            continue
    return "unknown"


def run(cmd: list[str], *, cwd: Path | str | None = None, timeout: int = 3600,
        log_path: Path | str | None = None, tool: str = "") -> subprocess.CompletedProcess:
    """Запустить внешнюю команду безопасно.

    Возвращает CompletedProcess при коде 0; иначе :class:`ToolRunError`.
    stdout/stderr при указании ``log_path`` пишутся в файл.
    """
    cmd = [str(c) for c in cmd]
    tool = tool or Path(cmd[0]).name
    log.debug("run: %s", " ".join(cmd))
    try:
        proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None,
                              capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as e:
        tail = (e.stderr or e.stdout or "")[-2000:] if isinstance(e.stderr, str) else ""
        raise ToolRunError(tool, cmd, None, tail) from e
    if log_path is not None:
        Path(log_path).write_text((proc.stdout or "") + (proc.stderr or ""), encoding="utf-8")
    if proc.returncode != 0:
        tail = ((proc.stdout or "") + (proc.stderr or ""))[-2000:]
        raise ToolRunError(tool, cmd, proc.returncode, tail)
    return proc
