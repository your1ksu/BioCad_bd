#!/usr/bin/env python3
from SHM_cold_spots_hot_spots import (
    find_motifs,
    HOTSPOT_MOTIFS,
    COLDSPOT_MOTIFS,
)

GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"


def run_test(test_name, test_func):
    """Вспомогательная функция для запуска и логирования тестов."""
    try:
        test_func()
        print(f"{GREEN}[УСПЕШНО]{RESET} {test_name}")
        return True
    except AssertionError as e:
        print(f"{RED}[ОШИБКА]{RESET} {test_name}")
        print(f"  -> Проверка не прошла (AssertionError)")
        return False
    except Exception as e:
        print(f"{RED}[КРАХ]{RESET} {test_name}")
        print(f"  -> Неожиданное исключение: {e}")
        return False


# --- СПИСОК ТЕСТОВ ---

def test_find_hotspots_clear():
    """1. Базовый поиск горячих точек (без гэпов)"""
    seq = "CGATAAATCG"
    matches = find_motifs(seq, HOTSPOT_MOTIFS, "hotspot")
    
    assert len(matches) == 1
    assert matches[0].motif_name == "TAA"
    assert matches[0].sequence == "TAA"
    assert matches[0].start_index == 3
    assert matches[0].end_index == 6
    assert matches[0].gapped is False


def test_find_coldspots_clear():
    """2. Базовый поиск холодных точек (без гэпов)"""
    seq = "ATGCCCGAT"
    matches = find_motifs(seq, COLDSPOT_MOTIFS, "coldspot")
    
    assert len(matches) == 2
    assert matches[0].motif_name == "SYC"
    assert matches[0].sequence == "GCC"
    assert matches[1].sequence == "CCC"


def test_find_motifs_with_gaps():
    """3. Поиск сквозь гэпы (проверка флага gapped)"""
    seq = "CGT-AATCG"
    matches = find_motifs(seq, HOTSPOT_MOTIFS, "hotspot")
    
    assert len(matches) == 1
    assert matches[0].motif_name == "TAA"
    assert matches[0].sequence == "TAA"
    assert matches[0].start_index == 2
    assert matches[0].end_index == 6
    assert matches[0].gapped is True


def test_overlapping_motifs():
    """4. Поиск перекрывающихся мотивов"""
    seq = "TAAATAA"
    matches = find_motifs(seq, HOTSPOT_MOTIFS, "hotspot")
    
    taa_matches = [m for m in matches if m.motif_name == "TAA"]
    assert len(taa_matches) == 2
    assert taa_matches[0].start_index == 0
    assert taa_matches[1].start_index == 4


def test_no_infinite_loop_on_trailing_gaps():
    """5. Защита от бесконечного цикла на хвостах из гэпов"""
    seq = "ATCG--------"
    matches = find_motifs(seq, HOTSPOT_MOTIFS, "hotspot")
    assert len(matches) == 0


# --- ГЛАВНЫЙ БЛОК ЗАПУСКА ---

if __name__ == "__main__":
    print("Запуск тестов скрипта SHM (без внешних библиотек)...\n")
    
    tests = [
        ("Поиск горячих точек (без гэпов)", test_find_hotspots_clear),
        ("Поиск холодных точек (без гэпов)", test_find_coldspots_clear),
        ("Поиск мотивов с гэпами", test_find_motifs_with_gaps),
        ("Поиск перекрывающихся мотивов", test_overlapping_motifs),
        ("Защита от бесконечного цикла", test_no_infinite_loop_on_trailing_gaps),
    ]
    
    passed_count = 0
    for name, func in tests:
        if run_test(name, func):
            passed_count += 1
            
    print(f"\nИтог: успешно пройдено {passed_count} из {len(tests)} тестов.")
    
    if passed_count != len(tests):
        exit(1)
    else:
        exit(0)
