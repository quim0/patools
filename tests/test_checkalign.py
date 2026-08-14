import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from checkalign.checkalign import (
    check_score_affine,
    check_score_affine2p,
    check_score_edit,
    parse_cigar,
    parse_sequence_pairs,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def run_checkalign(arguments, stdin=""):
    return subprocess.run(
        [sys.executable, "-m", "checkalign.checkalign", *arguments],
        cwd=REPOSITORY_ROOT,
        input=stdin,
        capture_output=True,
        text=True,
        check=False,
    )


class ParseCigarTests(unittest.TestCase):
    def test_parses_supported_operations(self):
        operations, repetitions = parse_cigar("12M2X3I4D")

        self.assertEqual(operations, ["M", "X", "I", "D"])
        self.assertEqual(repetitions, [12, 2, 3, 4])

    def test_rejects_incomplete_or_non_positive_cigars(self):
        for cigar in ("", "M", "0M", "-1M", "1M2", "1M1S", "1m"):
            with self.subTest(cigar=cigar):
                with self.assertRaises(ValueError):
                    parse_cigar(cigar)


class SequencePairTests(unittest.TestCase):
    def test_parses_sequence_markers(self):
        self.assertEqual(
            parse_sequence_pairs([">AC\n", "<AG\n"]),
            [("AC", "AG")],
        )

    def test_rejects_incomplete_or_unmarked_pairs(self):
        for lines in ([">AC\n"], ["AC\n", "AG\n"], ["<AC\n", ">AG\n"]):
            with self.subTest(lines=lines):
                with self.assertRaises(ValueError):
                    parse_sequence_pairs(lines)


class ScoringTests(unittest.TestCase):
    def test_edit_score(self):
        self.assertEqual(
            check_score_edit(5, ["M", "X", "I", "D"], [10, 2, 1, 2]),
            (True, 5),
        )

    def test_affine_score(self):
        self.assertEqual(
            check_score_affine(-7, ["M", "X", "I"], [2, 1, 2], 0, -1, -2, -2),
            (True, -7),
        )

    def test_affine2p_uses_the_better_gap_model(self):
        self.assertEqual(
            check_score_affine2p(
                -4, ["I"], [3], 0, -1, -1, -1, -10, -1
            ),
            (True, -4),
        )


class CheckalignCliTests(unittest.TestCase):
    def test_empty_input_is_a_clean_failure(self):
        result = run_checkalign(["-q", "-"])

        self.assertEqual(result.returncode, 1)
        self.assertNotIn("Traceback", result.stderr)

    def test_malformed_cigar_is_not_accepted(self):
        result = run_checkalign(["-q", "-"], "0 1M2\n")

        self.assertEqual(result.returncode, 1)
        self.assertIn("Invalid CIGAR at line 1", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_missing_cigar_count_is_a_clean_failure(self):
        result = run_checkalign(["-q", "-"], "0 M\n")

        self.assertEqual(result.returncode, 1)
        self.assertNotIn("Traceback", result.stderr)

    def test_three_field_row_is_rejected(self):
        result = run_checkalign(["-q", "-"], "0 1M A\n")

        self.assertEqual(result.returncode, 1)
        self.assertIn("Invalid result row at line 1", result.stderr)

    def test_embedded_sequence_failure_is_preserved(self):
        result = run_checkalign(["-q", "-"], "0 1M A T\n")

        self.assertEqual(result.returncode, 1)
        self.assertNotIn("Traceback", result.stderr)

    def test_quiet_valid_run_has_no_stdout(self):
        result = run_checkalign(["-q", "-"], "0 1M A A\n")

        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")

    def test_ground_truth_mismatch_is_reported_without_crashing(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            results_file = directory / "results.out"
            ground_truth_file = directory / "truth.out"
            results_file.write_text("0 1M\n")
            ground_truth_file.write_text("1 1X\n")

            result = run_checkalign(
                ["-q", "--ground-truth", str(ground_truth_file), str(results_file)]
            )

        self.assertEqual(result.returncode, 1)
        self.assertNotIn("Traceback", result.stderr)

    def test_sequence_pairs_are_reused_for_each_results_file(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            first_results = directory / "first.out"
            second_results = directory / "second.out"
            sequences = directory / "input.seq"
            first_results.write_text("0 1M\n")
            second_results.write_text("0 1M\n")
            sequences.write_text(">A\n<A\n")

            result = run_checkalign(
                [
                    "-q",
                    "--sequences",
                    str(sequences),
                    str(first_results),
                    str(second_results),
                ]
            )

        self.assertEqual(result.returncode, 0, result.stderr)

    def test_multiple_files_can_be_checked_in_parallel(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            first_results = directory / "first.out"
            second_results = directory / "second.out"
            first_results.write_text("0 1M A A\n")
            second_results.write_text("1 1X A T\n")

            result = run_checkalign(
                ["-q", "--jobs", "2", str(first_results), str(second_results)]
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")

    def test_jobs_must_be_positive(self):
        result = run_checkalign(["-q", "--jobs", "0", "-"])

        self.assertEqual(result.returncode, 2)
        self.assertNotIn("Traceback", result.stderr)

    def test_invalid_penalties_produce_an_argument_error(self):
        result = run_checkalign(["-q", "--penalties", "bad", "-"])

        self.assertEqual(result.returncode, 2)
        self.assertNotIn("Traceback", result.stderr)


if __name__ == "__main__":
    unittest.main()
