"""
Тесты для scripts/tree_analytics/tree_analytics.py

Проверяем только extract_statistics — чистую функцию сбора статистики.
plot_results не тестируем (она рисует картинку, а не возвращает данные).
"""
from tree_analytics import extract_statistics


def _sample_data():
    return {
        "geneA": {
            "mrbayes": {
                "clades": [
                    {
                        "size": 3,
                        "depth": 1,
                        "ancestor_to_leaves": {"max": 0.1, "mean": 0.05, "min": 0.01},
                    },
                    {
                        "size": 5,
                        "depth": 2,
                        "ancestor_to_leaves": {"max": 0.2, "mean": 0.10, "min": 0.02},
                    },
                ]
            }
        },
        "geneB": {
            "mrbayes": {
                "clades": [
                    {
                        "size": 2,
                        "depth": 0,
                        "ancestor_to_leaves": {"max": 0.05, "mean": 0.03, "min": 0.01},
                    },
                ]
            }
        },
    }


def test_extract_statistics_genes_and_counts():
    genes, clades_counts, sizes, depths, ancestor_distances = extract_statistics(_sample_data())

    assert genes == ["geneA", "geneB"]
    assert clades_counts == [2, 1]


def test_extract_statistics_sizes():
    _, _, sizes, _, _ = extract_statistics(_sample_data())
    assert sizes == [3, 5, 2]


def test_extract_statistics_depths_from_field():
    # depth теперь берётся напрямую из готового поля, а не вычисляется
    _, _, _, depths, _ = extract_statistics(_sample_data())
    assert depths == [1, 2, 0]


def test_extract_statistics_ancestor_distances():
    _, _, _, _, ancestor_distances = extract_statistics(_sample_data())
    assert ancestor_distances == [0.05, 0.10, 0.03]


def test_extract_statistics_missing_depth_defaults_to_zero():
    data = {
        "geneC": {
            "mrbayes": {
                "clades": [{"size": 4}]  # нет ни depth, ни ancestor_to_leaves
            }
        }
    }
    genes, clades_counts, sizes, depths, ancestor_distances = extract_statistics(data)

    assert sizes == [4]
    assert depths == [0]
    assert ancestor_distances == []  # нет данных - не добавляем, а не подставляем 0


def test_extract_statistics_no_clades():
    data = {"geneD": {"mrbayes": {"clades": []}}}
    genes, clades_counts, sizes, depths, ancestor_distances = extract_statistics(data)

    assert genes == ["geneD"]
    assert clades_counts == [0]
    assert sizes == []
    assert depths == []
    assert ancestor_distances == []


def test_extract_statistics_empty_data():
    genes, clades_counts, sizes, depths, ancestor_distances = extract_statistics({})
    assert genes == []
    assert clades_counts == []
    assert sizes == []
    assert depths == []
    assert ancestor_distances == []
