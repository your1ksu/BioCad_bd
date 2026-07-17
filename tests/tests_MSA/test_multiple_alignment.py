from pathlib import Path
from tempfile import TemporaryDirectory
import sys
import unittest
from unittest.mock import patch


MSA_DIR = Path(__file__).resolve().parents[2] / "scripts" / "03_multiple_alignment"
sys.path.insert(0, str(MSA_DIR))

import multiple_alignment


class TestMsaFinal(unittest.TestCase):
    def test_discover_fasta_files_recursively(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            tmp_path = Path(temporary_directory)
            input_dir = tmp_path / "input"
            nested_dir = input_dir / "vj"
            nested_dir.mkdir(parents=True)

            fasta_file = input_dir / "root.fasta"
            nested_fasta_file = nested_dir / "nested.fa"
            ignored_file = input_dir / "notes.txt"

            fasta_file.write_text(">seq1\nACGT\n", encoding="utf-8")
            nested_fasta_file.write_text(">seq2\nAGT\n", encoding="utf-8")
            ignored_file.write_text("not fasta\n", encoding="utf-8")

            result = multiple_alignment.discover_fasta_files(input_dir)

            self.assertEqual(result, sorted([fasta_file, nested_fasta_file]))

    def test_find_mafft_uses_explicit_path(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            mafft_bin = Path(temporary_directory) / "mafft"
            mafft_bin.write_text("#!/bin/sh\n", encoding="utf-8")

            self.assertEqual(multiple_alignment.find_mafft(mafft_bin), mafft_bin)

    def test_run_mafft_writes_alignment(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            tmp_path = Path(temporary_directory)
            input_fasta = tmp_path / "sequences.fasta"
            output_fasta = tmp_path / "sequences_aligned.fasta"
            input_fasta.write_text(">seq1\nACGT\n>seq2\nAGT\n", encoding="utf-8")

            aligned_fasta = ">seq1\nACGT\n>seq2\nA-GT\n"

            def fake_mafft(command, stdout, check):
                self.assertEqual(
                    command,
                    ["mafft", "--auto", "--quiet", "--preservecase", str(input_fasta)],
                )
                self.assertTrue(check)
                stdout.write(aligned_fasta)

            with patch("multiple_alignment.subprocess.run", side_effect=fake_mafft):
                multiple_alignment.run_mafft(Path("mafft"), input_fasta, output_fasta)

            self.assertEqual(output_fasta.read_text(encoding="utf-8"), aligned_fasta)

    def test_main_aligns_files_and_writes_manifest(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            tmp_path = Path(temporary_directory)
            input_dir = tmp_path / "input"
            output_dir = tmp_path / "output"
            nested_dir = input_dir / "v"
            nested_dir.mkdir(parents=True)

            input_fasta = nested_dir / "group1.fasta"
            input_fasta.write_text(">seq1\nACGT\n>seq2\nAGT\n", encoding="utf-8")

            mafft_bin = tmp_path / "mafft"
            mafft_bin.write_text("#!/bin/sh\n", encoding="utf-8")

            test_argv = [
                "multiple_alignment.py",
                "--input",
                str(input_dir),
                "--output",
                str(output_dir),
                "--mafft",
                str(mafft_bin),
            ]

            def fake_mafft(command, stdout, check):
                self.assertEqual(command[0], str(mafft_bin))
                self.assertEqual(command[-1], str(input_fasta))
                self.assertTrue(check)
                stdout.write(">seq1\nACGT\n>seq2\nA-GT\n")

            with patch.object(sys, "argv", test_argv), patch(
                "multiple_alignment.subprocess.run", side_effect=fake_mafft
            ):
                multiple_alignment.main()

            aligned_fasta = output_dir / "v" / "group1_aligned.fasta"
            manifest = output_dir / "manifest.tsv"

            self.assertEqual(
                aligned_fasta.read_text(encoding="utf-8"),
                ">seq1\nACGT\n>seq2\nA-GT\n",
            )
            self.assertEqual(
                manifest.read_text(encoding="utf-8"),
                "input_fasta\taligned_fasta\nv/group1.fasta\tv/group1_aligned.fasta\n",
            )


if __name__ == "__main__":
    unittest.main()
