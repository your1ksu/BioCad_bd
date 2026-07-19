"""
Тесты для scripts/translate_to_amino/translate_to_amino.py

Запуск:  pytest tests/test_translate_to_amino.py -v
"""
import os
import sys

import pytest

SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS_DIR)

import translate_to_amino as tta  # noqa: E402


# ---------- trim_to_complete_codons ----------

def test_trim_to_complete_codons_exact_multiple_of_three():
    assert tta.trim_to_complete_codons("ATGAAATTT", 0) == "ATGAAATTT"


def test_trim_to_complete_codons_drops_trailing_leftover():
    # 11 нуклеотидов, рамка 0 -> 9 целых + 2 "хвостовых" отбрасываются
    assert tta.trim_to_complete_codons("ATGAAATTTGG", 0) == "ATGAAATTT"


def test_trim_to_complete_codons_applies_frame_offset():
    # сдвиг на 1 нуклеотид вправо перед тем, как резать по три
    assert tta.trim_to_complete_codons("XATGAAATTT", 1) == "ATGAAATTT"


# ---------- translate_frame ----------

def test_translate_frame_basic():
    # ATG AAA TAA -> M K *
    assert tta.translate_frame("ATGAAATAA", 0) == "MK*"


def test_translate_frame_too_short_returns_none():
    assert tta.translate_frame("AT", 0) is None


def test_translate_frame_with_offset():
    # сдвиг на 1: "TGAAATAA" -> TGA AAT AA(отброшено) -> "*N"
    assert tta.translate_frame("XTGAAATAA", 1) == "*N"


# ---------- classify_stop_codon ----------

def test_classify_stop_codon_no_stop():
    assert tta.classify_stop_codon("MKTA") == "no_stop"


def test_classify_stop_codon_stop_at_end():
    assert tta.classify_stop_codon("MKTA*") == "stop_at_end"


def test_classify_stop_codon_premature():
    assert tta.classify_stop_codon("MK*TA") == "premature_stop"


def test_classify_stop_codon_multiple_stops():
    assert tta.classify_stop_codon("MK*TA*") == "premature_stop"


# ---------- strip_trailing_stop ----------

def test_strip_trailing_stop_removes_star():
    assert tta.strip_trailing_stop("MKTA*") == "MKTA"


def test_strip_trailing_stop_no_star_unchanged():
    assert tta.strip_trailing_stop("MKTA") == "MKTA"


# ---------- read_aa_reference ----------

def test_read_aa_reference_basic(tmp_path):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">acc|IGHV1-2*01|extra\nMKTAYIRQKS\n>acc2|IGHJ1*01|extra\nWGQGTLVTVSS\n")
    ref = tta.read_aa_reference(str(fasta))
    assert ref == {"IGHV1-2*01": "MKTAYIRQKS", "IGHJ1*01": "WGQGTLVTVSS"}


# ---------- best_reference_match / is_confident_match ----------

def test_best_reference_match_finds_exact():
    reference = {"GENE_A": "MKTAYIRQKS", "GENE_B": "WGQGTLVTVSS"}
    name, score, ref_len = tta.best_reference_match("MKTAYIRQKS", reference)
    assert name == "GENE_A"
    assert ref_len == 10
    assert score == 10  # полное совпадение


def test_is_confident_match_uses_reference_length_not_protein_length():
    # белок длиннее гена из справочника, но полностью его содержит -> должно засчитаться
    reference = {"GENE_A": "MKTAYIRQKS"}
    protein = "MKTAYIRQKS" + "WGQGTLVTVSS"  # V-ген + J-ген склеены
    _, score, ref_len = tta.best_reference_match(protein, reference)
    assert tta.is_confident_match(score, ref_len) is True


def test_is_confident_match_rejects_unrelated_sequence():
    reference = {"GENE_A": "MKTAYIRQKS"}
    _, score, ref_len = tta.best_reference_match("CCCCCCCCCC", reference)
    assert tta.is_confident_match(score, ref_len) is False


# ---------- process_sequence (сквозная логика выбора рамки + сверки) ----------

CODON = {
    "M": "ATG", "K": "AAA", "T": "ACC", "A": "GCC", "Y": "TAC", "I": "ATC",
    "R": "CGC", "Q": "CAG", "S": "TCC", "W": "TGG", "G": "GGC", "L": "CTC",
    "V": "GTC",
}
STOP = "TAA"


def _build_nt(protein, trailing_stop=True, prefix=""):
    nt = prefix + "".join(CODON[a] for a in protein)
    if trailing_stop:
        nt += STOP
    return nt


@pytest.fixture
def reference():
    return {
        "IGHV-TEST*01": "MKTAYIRQKS",
        "IGHJ-TEST*01": "WGQGTLVTVSS",
    }


def test_process_sequence_frame_zero_works(reference):
    full_protein = "MKTAYIRQKSWGQGTLVTVSS"
    nt = _build_nt(full_protein)
    counters = {"no_stop": 0, "stop_at_end": 0, "matched": 0, "no_match": 0, "premature_stop": 0}

    status, protein, frame = tta.process_sequence(nt, reference, counters)

    assert status == "ok"
    assert protein == full_protein
    assert frame == 1
    assert counters["matched"] == 1
    assert counters["stop_at_end"] == 1


def test_process_sequence_needs_frame_shift(reference):
    full_protein = "MKTAYIRQKSWGQGTLVTVSS"
    nt = "A" + _build_nt(full_protein)  # сдвигаем рамку на 1 лишним нуклеотидом спереди
    counters = {"no_stop": 0, "stop_at_end": 0, "matched": 0, "no_match": 0, "premature_stop": 0}

    status, protein, frame = tta.process_sequence(nt, reference, counters)

    assert status == "ok"
    assert protein == full_protein
    assert frame == 2  # потребовалась рамка 2


def test_process_sequence_no_match_when_unrelated_to_reference(reference):
    nt = _build_nt("MKTAYIRQKS")  # валидно по стоп-кодону, но это не то, что в reference
    counters = {"no_stop": 0, "stop_at_end": 0, "matched": 0, "no_match": 0, "premature_stop": 0}

    # намеренно берём последовательность, которая совпадёт (это часть V-гена) —
    # проверим обратный случай: заведомо непохожая последовательность
    nt_unrelated = _build_nt("GGGGGGGGGG")
    status, protein, frame = tta.process_sequence(nt_unrelated, reference, counters)

    assert status == "no_match"
    assert counters["no_match"] == 1


def test_process_sequence_premature_stop_in_all_frames(reference):
    # стоп-кодон специально не в конце, и другие рамки тоже не спасают
    nt = "ATG" + STOP + "AAAGGCTTT"  # M * K G F -> стоп на 2-й позиции
    counters = {"no_stop": 0, "stop_at_end": 0, "matched": 0, "no_match": 0, "premature_stop": 0}

    status, protein, frame = tta.process_sequence(nt, reference, counters)

    # либо ни одна рамка не даёт валидный стоп (premature_stop),
    # либо валидная рамка есть, но не совпадает со справочником (no_match) —
    # в обоих случаях protein должен быть None
    assert status in ("premature_stop", "no_match")
    assert protein is None


def test_process_sequence_without_reference_accepts_first_valid_frame():
    full_protein = "MKTAYIRQKS"
    nt = _build_nt(full_protein)
    counters = {"no_stop": 0, "stop_at_end": 0, "matched": 0, "no_match": 0, "premature_stop": 0}

    status, protein, frame = tta.process_sequence(nt, {}, counters)

    assert status == "ok"
    assert protein == full_protein
    assert counters["matched"] == 1


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
