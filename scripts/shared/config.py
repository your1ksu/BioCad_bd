from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class PipelineConfig:
    # ---------- paths ----------
    ref_dir: Path = Path("data/references")          # папка с germline референсами (IGHV.fasta, IGHD.fasta, ...)
    report_dir: Optional[Path] = None                # куда складывать результаты (авто-генерация если None)

    # ---------- filter (step 1) ----------
    min_junction: int = 15                           # мин. длина junction-региона для фильтрации
    max_junction: int = 300                          # макс. длина junction
    v_min_fraction: float = 0.5                      # мин. доля V-гена в последовательности
    filter_mode: str = "germline"                    # режим фильтрации ("germline" | ...)

    # ---------- group_by_germlines (step 2) ----------
    grouping_strategy: str = "gene"                  # allele (IGHV1-2*01_IGHJ4*01), gene (IGHV1-2_IGHJ4), v_only (IGHV1-2)

    # ---------- filter_by_symbol_count (step 2b) ----------
    min_group_size: int = 5                          # мин. число последовательностей в группе (меньшие удаляются)
    max_group_size: int = 100                        # макс. число последовательностей в группе (большие обрезаются)

    # ---------- build_trees_iqtree (step 4a) ----------
    iqtree_model: str = "GTR+F+I+G4"                # модель нуклеотидных замен для IQ-TREE (GTR+F+I+G4 — быстрее MFP)

    # ---------- build_trees_mrbayes (step 4b) ----------
    mb_ngen_default: int = 200_000                   # число генераций MCMC (если не рассчитывается динамически)
    mb_ngen_min: int = 50_000                        # минимум генераций при динамическом расчёте
    mb_ngen_max: int = 500_000                       # максимум генераций при динамическом расчёте
    mb_burnin_frac: float = 0.25                     # доля burn-in (первые итерации, которые выкидываются)

    # ---------- clade_search / tree_analytics ----------
    posterior_min: float = 0.95                      # мин. posterior probability для confident clade
    ufboot_min: float = 95.0                         # мин. UFBoot support (%) для confident clade
    alrt_min: float = 80.0                           # мин. SH-aLRT support (%) для confident clade

    # ---------- parallel execution ----------
    parallel_trees: bool = False                     # запускать IQ-TREE и MrBayes параллельно (в 2 потока)
    workers: int = 0                                 # число параллельных процессов (0 = все доступные ядра)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2, default=str)

    def to_dict(self) -> dict:
        d = {}
        for k, v in self.__dataclass_fields__.items():
            val = getattr(self, k)
            if isinstance(val, Path):
                val = str(val)
            d[k] = val
        return d

    @classmethod
    def load(cls, path: Path) -> PipelineConfig:
        with open(path) as f:
            data = json.load(f)
        cfg = cls()
        for k, v in data.items():
            if hasattr(cfg, k):
                if k in ("ref_dir", "report_dir") and v:
                    v = Path(v)
                setattr(cfg, k, v)
        return cfg


DEFAULT_CONFIG = PipelineConfig()
