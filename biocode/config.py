"""Конфигурация прогона (единый источник настроек).

``RunConfig`` собирает все параметры пайплайна: выбор рабочей последовательности,
правила группировки, режим MAFFT, настройки IQ-TREE/MrBayes, параллелизм, сиды,
пути к инструментам и веса метрик ранжирования. Сериализуется в манифест прогона
для воспроизводимости.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path

from .errors import ConfigError

# какую колонку AIRR брать как рабочую последовательность
WORKING_SEQ_CHOICES = ("sequence_alignment", "sequence_vdj", "sequence")
GROUP_BY_CHOICES = ("v_j_gene", "v_gene", "v_j_call", "v_call")


@dataclass
class RunConfig:
    # ── вход/выход ──────────────────────────────────────────────────────────
    input: str = ""                         # путь к AIRR-TSV или FASTA
    out: str = "EDU/results/biocode_run"    # каталог прогона
    working_seq: str = "sequence_alignment"  # см. WORKING_SEQ_CHOICES
    # ── группировка ─────────────────────────────────────────────────────────
    group_by: str = "v_j_gene"              # см. GROUP_BY_CHOICES
    min_group_size: int = 4                 # меньше — дерево неинформативно (skip)
    locus: str = "IGH"                      # фильтр локуса (пусто = не фильтровать)
    productive_only: bool = True
    limit_groups: int = 0                   # 0 = без лимита (для быстрых прогонов)
    # ── выравнивание ────────────────────────────────────────────────────────
    mafft_mode: str = "linsi"               # linsi | auto | ginsi | einsi
    # ── ML (IQ-TREE) ────────────────────────────────────────────────────────
    iqtree_model: str = "MFP"               # ModelFinder Plus
    ufboot: int = 1000                      # -B
    alrt: int = 1000                        # -alrt
    asr: bool = True                        # --asr (реконструкция предков)
    # ── Bayes (MrBayes) ─────────────────────────────────────────────────────
    run_bayes: bool = False                 # включать байесовский путь
    mb_ngen: int = 200000                   # длина цепи MCMC
    mb_burnin_frac: float = 0.25
    # ── ресурсы/детерминизм ─────────────────────────────────────────────────
    make_plots: bool = True                 # рисовать дерево + FR/CDR-полосу на группу
    threads: str = "AUTO"                   # -T для IQ-TREE
    jobs: int = 1                           # групп в параллель
    seed: int = 12345
    timeout_s: int = 3600                   # таймаут на внешний вызов
    resume: bool = True                     # пропускать уже посчитанные группы
    # ── явные пути к бинарям (пусто = автопоиск) ────────────────────────────
    mafft_bin: str = ""
    iqtree_bin: str = ""
    mrbayes_bin: str = ""
    igblast_bin: str = ""
    # ── germline-референс (для outgroup/фолбэка) ────────────────────────────
    imgt_ref_dir: str = "structural_platform/IMGT_V-QUEST_reference_directory"
    # ── веса метрик ранжирования (блок B9) ──────────────────────────────────
    weights: dict[str, float] = field(default_factory=lambda: {
        "recurrence": 2.0, "region_cdr": 1.5, "replacement": 1.0,
        "support": 1.0, "persistence": 1.0,
    })
    log_level: str = "INFO"

    def validate(self) -> "RunConfig":
        if self.working_seq not in WORKING_SEQ_CHOICES:
            raise ConfigError(f"working_seq='{self.working_seq}', допустимо: {WORKING_SEQ_CHOICES}")
        if self.group_by not in GROUP_BY_CHOICES:
            raise ConfigError(f"group_by='{self.group_by}', допустимо: {GROUP_BY_CHOICES}")
        if self.min_group_size < 2:
            raise ConfigError("min_group_size должен быть ≥ 2")
        if not (0.0 <= self.mb_burnin_frac < 1.0):
            raise ConfigError("mb_burnin_frac должен быть в [0,1)")
        if self.jobs < 1:
            raise ConfigError("jobs должен быть ≥ 1")
        return self

    @property
    def out_path(self) -> Path:
        return Path(self.out)

    @property
    def run_dir(self) -> Path:
        return self.out_path / "run"

    def to_dict(self) -> dict:
        return asdict(self)
