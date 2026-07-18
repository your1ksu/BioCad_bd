"""BIOCODE core — движок филогенетики и анализа мутаций антител.

Публичный API (стабильные точки входа). Подробности — в
``EDU/docs/PHYLO_ARCHITECTURE.md`` и ``EDU/docs/PHYLO_PLAN.md``.
"""
from __future__ import annotations

__version__ = "0.1.0"

from .errors import (
    BiocodeError, ConfigError, InputError, ValidationError,
    ToolNotFoundError, ToolRunError,
)
from .config import RunConfig
from .model import SequenceRecord, Group, TreeResult, Mutation

__all__ = [
    "__version__",
    "BiocodeError", "ConfigError", "InputError", "ValidationError",
    "ToolNotFoundError", "ToolRunError",
    "RunConfig", "SequenceRecord", "Group", "TreeResult", "Mutation",
]
