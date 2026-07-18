"""Структурное логирование прогона.

Единая точка получения логгера. По умолчанию пишет в консоль; при указании
файла — дублирует в файл прогона (`<out>/run/run.log`). Формат единообразный,
с временными метками, чтобы логи групп были сопоставимы.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

_CONFIGURED = False
_FMT = "%(asctime)s %(levelname)-7s %(name)s | %(message)s"
_DATEFMT = "%H:%M:%S"


def configure(level: str = "INFO", logfile: Path | str | None = None) -> None:
    """Однократная настройка корневого логгера пакета."""
    global _CONFIGURED
    root = logging.getLogger("biocode")
    root.setLevel(getattr(logging, level.upper(), logging.INFO))
    # избегаем дублирования хендлеров при повторных вызовах
    for h in list(root.handlers):
        root.removeHandler(h)
    fmt = logging.Formatter(_FMT, _DATEFMT)
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    root.addHandler(sh)
    if logfile is not None:
        Path(logfile).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(logfile, encoding="utf-8")
        fh.setFormatter(fmt)
        root.addHandler(fh)
    root.propagate = False
    _CONFIGURED = True


def get_logger(name: str = "biocode") -> logging.Logger:
    """Логгер пакета. Настраивает дефолт, если configure() ещё не звали."""
    if not _CONFIGURED:
        configure()
    full = name if name.startswith("biocode") else f"biocode.{name}"
    return logging.getLogger(full)
