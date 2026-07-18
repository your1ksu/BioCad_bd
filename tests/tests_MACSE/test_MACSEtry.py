from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import sys
import unittest
from unittest.mock import patch


MACSE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MACSE_DIR))

import MACSEtry


class TestMacseTry(unittest.TestCase):
    def test_discover_fasta_files_recursively(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            tmp_path = Path(temporary_directory)
            input_dir = tmp_path / "input"
            nested_dir = input_dir / "v"
            nested_dir.mkdir(parents=True)

            fasta_file = input_dir / "root.fasta"
            nested_fasta_file = nested_dir / "nested.fa"
            uppercase_fasta_file = nested_dir / "upper.FNA"
            ignored_file = input_dir / "notes.txt"

            fasta_file.write_text(">seq1\nACGT\n", encoding="utf-8")
            nested_fasta_file.write_text(">seq2\nAGT\n", encoding="utf-8")
            uppercase_fasta_file.write_text(">seq3\nATG\n", encoding="utf-8")
            ignored_file.write_text("not fasta\n", encoding="utf-8")

            result = MACSEtry.discover_fasta_files(input_dir)

            self.assertEqual(
                result,
                sorted([fasta_file, nested_fasta_file, uppercase_fasta_file]),
            )

    def test_amino_acid_output_path_adds_aa_suffix(self) -> None:
        output_fasta = Path("group1_aligned.fasta")

        self.assertEqual(
            MACSEtry.amino_acid_output_path(output_fasta),
            Path("group1_aligned_aa.fasta"),
        )

    def test_find_macse_uses_explicit_path(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            macse_bin = Path(temporary_directory) / "macse"
            macse_bin.write_text("#!/bin/sh\n", encoding="utf-8")

            self.assertEqual(MACSEtry.find_macse(macse_bin), macse_bin)

    def test_find_macse_uses_path_lookup(self) -> None:
        with patch("MACSEtry.shutil.which", return_value="/usr/bin/macse"):
            self.assertEqual(MACSEtry.find_macse(None), Path("/usr/bin/macse"))

    def test_find_macse_raises_when_not_found(self) -> None:
        with patch("MACSEtry.shutil.which", return_value=None):
            with self.assertRaises(FileNotFoundError):
                MACSEtry.find_macse(None)

    def test_run_macse_writes_nt_and_aa_alignments(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            tmp_path = Path(temporary_directory)
            input_fasta = tmp_path / "sequences.fasta"
            output_fasta = tmp_path / "sequences_aligned.fasta"
            aa_output_fasta = tmp_path / "sequences_aligned_aa.fasta"
            input_fasta.write_text(">seq1\nATG\n>seq2\nATG\n", encoding="utf-8")

            def fake_macse(command, check):
                self.assertEqual(command[0], "macse")
                self.assertEqual(command[1:4], ["-prog", "alignSequences", "-seq"])
                self.assertEqual(command[4], str(input_fasta))
                self.assertEqual(command[5], "-out_NT")
                self.assertEqual(command[7], "-out_AA")
                self.assertTrue(check)

                Path(command[6]).write_text(">seq1\nATG\n>seq2\nATG\n", encoding="utf-8")
                Path(command[8]).write_text(">seq1\nM\n>seq2\nM\n", encoding="utf-8")

            with patch("MACSEtry.subprocess.run", side_effect=fake_macse):
                MACSEtry.run_macse(Path("macse"), input_fasta, output_fasta)

            self.assertEqual(
                output_fasta.read_text(encoding="utf-8"),
                ">seq1\nATG\n>seq2\nATG\n",
            )
            self.assertEqual(
                aa_output_fasta.read_text(encoding="utf-8"),
                ">seq1\nM\n>seq2\nM\n",
            )
            self.assertFalse((tmp_path / "sequences_aligned.tmp.fasta").exists())
            self.assertFalse((tmp_path / "sequences_aligned_aa.tmp.fasta").exists())

    def test_run_macse_cleans_partial_files_after_error(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            tmp_path = Path(temporary_directory)
            input_fasta = tmp_path / "sequences.fasta"
            output_fasta = tmp_path / "sequences_aligned.fasta"
            aa_output_fasta = tmp_path / "sequences_aligned_aa.fasta"
            input_fasta.write_text(">seq1\nATG\n", encoding="utf-8")
            output_fasta.write_text("partial nt\n", encoding="utf-8")
            aa_output_fasta.write_text("partial aa\n", encoding="utf-8")

            def fake_macse(command, check):
                Path(command[6]).write_text("temporary nt\n", encoding="utf-8")
                Path(command[8]).write_text("temporary aa\n", encoding="utf-8")
                raise subprocess.CalledProcessError(1, command)

            with patch("MACSEtry.subprocess.run", side_effect=fake_macse):
                with self.assertRaises(subprocess.CalledProcessError):
                    MACSEtry.run_macse(Path("macse"), input_fasta, output_fasta)

            self.assertFalse(output_fasta.exists())
            self.assertFalse(aa_output_fasta.exists())
            self.assertFalse((tmp_path / "sequences_aligned.tmp.fasta").exists())
            self.assertFalse((tmp_path / "sequences_aligned_aa.tmp.fasta").exists())

    def test_main_aligns_files_and_writes_manifest(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            tmp_path = Path(temporary_directory)
            input_dir = tmp_path / "input"
            output_dir = tmp_path / "output"
            nested_dir = input_dir / "v"
            nested_dir.mkdir(parents=True)

            input_fasta = nested_dir / "group1.fasta"
            input_fasta.write_text(">seq1\nATG\n>seq2\nATG\n", encoding="utf-8")

            macse_bin = tmp_path / "macse"
            macse_bin.write_text("#!/bin/sh\n", encoding="utf-8")

            test_argv = [
                "MACSEtry.py",
                "--input",
                str(input_dir),
                "--output",
                str(output_dir),
                "--macse",
                str(macse_bin),
            ]

            def fake_macse(command, check):
                self.assertEqual(command[0], str(macse_bin))
                self.assertEqual(command[4], str(input_fasta))
                self.assertTrue(check)
                Path(command[6]).write_text(">seq1\nATG\n>seq2\nATG\n", encoding="utf-8")
                Path(command[8]).write_text(">seq1\nM\n>seq2\nM\n", encoding="utf-8")

            with patch.object(sys, "argv", test_argv), patch(
                "MACSEtry.subprocess.run", side_effect=fake_macse
            ):
                MACSEtry.main()

            aligned_fasta = output_dir / "v" / "group1_aligned.fasta"
            aligned_aa_fasta = output_dir / "v" / "group1_aligned_aa.fasta"
            manifest = output_dir / "manifest.tsv"

            self.assertEqual(
                aligned_fasta.read_text(encoding="utf-8"),
                ">seq1\nATG\n>seq2\nATG\n",
            )
            self.assertEqual(
                aligned_aa_fasta.read_text(encoding="utf-8"),
                ">seq1\nM\n>seq2\nM\n",
            )
            self.assertEqual(
                manifest.read_text(encoding="utf-8"),
                (
                    "input_fasta\taligned_nt_fasta\taligned_aa_fasta\n"
                    "v/group1.fasta\tv/group1_aligned.fasta\t"
                    "v/group1_aligned_aa.fasta\n"
                ),
            )


if __name__ == "__main__":
    unittest.main()
