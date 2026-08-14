#!/usr/bin/env python3
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
import multiprocessing
import os
import re
import sys
from tqdm import tqdm
from rich.console import Console
from rich.table import Table

PROGRESS_CHECK = 5000
CIGAR_PATTERN = re.compile(r'(?:[1-9]\d*[MXID])+\Z')
CIGAR_ELEMENT_PATTERN = re.compile(r'([1-9]\d*)([MXID])')

console = Console()
error_console = Console(stderr=True)
_WORKER_CONFIG = None


@dataclass(frozen=True)
class CheckConfig:
    distance_function: str
    penalties: tuple
    with_mismatches: bool
    ground_truth_scores: tuple = None
    sequence_pairs: tuple = None
    collect_plot_scores: bool = False


@dataclass
class FileCheckResult:
    filename: str
    correct: int = 0
    incorrect: int = 0
    score_total: int = 0
    score_count: int = 0
    max_score: int = 0
    return_code: int = 0
    errors: list = field(default_factory=list)
    incorrect_rows: list = field(default_factory=list)
    plot_scores: list = field(default_factory=list)


def positive_integer(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return value

def print_report(correct, incorrect, results_file):
    console.print(f"([bold]{results_file}[/bold]) Correct=[green]{correct}[/green], Incorrect=[red]{incorrect}[/red], Accuracy=[bold]{correct*100/(correct+incorrect):.2f}%[/bold]")

def update_incorrect_cigars_table(table, line_num, score, cigar, cigar_score, gt_score=None):
    args = [str(line_num), str(score), cigar, str(cigar_score)]
    if len(table.columns) == 5:
        args.append(str(gt_score) if gt_score is not None else "N/A")
    table.add_row(*args)

def generate_incorrect_cigars_table(filename, with_gt):
    table = Table(title=f"Incorrect Alignments ({filename})")
    table.add_column("Alig. #")
    table.add_column("Score")
    table.add_column("CIGAR", max_width=80)
    table.add_column("CIGAR score")
    if with_gt:
        table.add_column("Ground truth score")

    return table


def parse_cigar(cigar):
    """Parse the supported CIGAR subset, rejecting partial or empty input."""
    if not CIGAR_PATTERN.fullmatch(cigar):
        raise ValueError(
            "CIGAR must contain positive counts followed by M, X, I, or D"
        )

    elements = CIGAR_ELEMENT_PATTERN.findall(cigar)
    cigar_reps = [int(reps) for reps, _ in elements]
    cigar_ops = [op for _, op in elements]
    return cigar_ops, cigar_reps


def parse_sequence_pairs(lines):
    """Parse pairs in the .seq format: a >pattern line and a <text line."""
    if len(lines) % 2 != 0:
        raise ValueError("Sequences file must contain complete pairs of lines")

    pairs = []
    for index in range(0, len(lines), 2):
        pattern_line = lines[index].strip()
        text_line = lines[index + 1].strip()
        pair_number = index // 2 + 1
        if not pattern_line.startswith('>') or not text_line.startswith('<'):
            raise ValueError(
                f"Invalid sequence pair {pair_number}: expected >pattern and <text"
            )
        pairs.append((pattern_line[1:], text_line[1:]))
    return pairs


def plot_cummulative_scores(data):
    import matplotlib.pyplot as plt
    import numpy as np

    try:
        plt.style.use('seaborn-v0_8-paper')
    except:
        pass

    plt.rcParams['figure.dpi']= 300

    fig, ax = plt.subplots()

    min_x = float('inf')
    max_x = 0

    for k, v in data.items():
        ydict = {}
        for score in v:
            if score not in ydict.keys():
                ydict[score] = 0
            ydict[score] += 1

        x = []
        y = []
        for idx, skey in enumerate(sorted(ydict.keys())):

            if skey < min_x:
                min_x = skey
            if skey > max_x:
                max_x = skey

            if True:
                x.append(skey)
                y.append(ydict[skey])
                if idx > 0:
                    y[idx] += y[idx-1]

        print(f"Plotting {k} with {len(x)} points")
        if k == 'ground_truth':
            label = 'Optimal solution'
        else:
            #label = k.split('/')[-2]
            # TODO: file name
            label = k
        ax.plot(x, y, label=label, linewidth=0.3)


    # add vertical line at x = 100, 500, 1000, 5000
    ax.axvline(x=100, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(x=500, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(x=1000, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(x=5000, color='gray', linestyle='--', linewidth=0.5)

    ax.set(xlabel='Score', ylabel='Cummulative count',
           title='Cummulative scores')
    #ax.set_xlim(min_x-1, max_x+1)
    ax.legend()
    plt.savefig('cummulative_scores.svg')

def check_score_edit(score, cigar_ops, cigar_reps):
    score_calc = 0
    for idx, op in enumerate(cigar_ops):
        reps = cigar_reps[idx]
        if op == 'M':
            continue
        elif op == 'X':
            score_calc += reps
        elif op in ['I', 'D']:
            score_calc += reps

    return (score == score_calc, score_calc)

def check_score_affine(score, cigar_ops, cigar_reps, M, X, O, E):
    score_calc = 0
    for idx, op in enumerate(cigar_ops):
        reps = cigar_reps[idx]
        if op == 'M':
            score_calc += M * reps
        elif op == 'X':
            score_calc += X * reps
        elif op in ['I', 'D']:
            score_calc += O + E * reps

    return (abs(score) == abs(score_calc), score_calc)

def check_score_affine2p(score, cigar_ops, cigar_reps, M, X, O1, E1, O2, E2):
    score_calc = 0
    for idx, op in enumerate(cigar_ops):
        reps = cigar_reps[idx]
        if op == 'M':
            score_calc += M * reps
        elif op == 'X':
            score_calc += X * reps
        elif op in ['I', 'D']:
            if O1 < 0 and E1 < 0 and O2 < 0 and E2 < 0:
                score_calc += max(
                    O1 + E1 * reps,
                    O2 + E2 * reps
                )
            else:
                score_calc += min(
                    O1 + E1 * reps,
                    O2 + E2 * reps
                )

    return (abs(score) == abs(score_calc), score_calc)

def check_cigar_sequences(score, cigar_ops, cigar_reps, pattern, text, with_mismatches=True):
    text_pos = 0
    pattern_pos = 0

    try:
        for idx, op in enumerate(cigar_ops):
            reps = cigar_reps[idx]
            for _ in range(reps):
                if op == 'M':
                    if with_mismatches and (pattern[pattern_pos] != text[text_pos]):
                        return False
                    pattern_pos += 1
                    text_pos += 1
                elif op == 'X':
                    if pattern[pattern_pos] == text[text_pos]:
                        return False
                    pattern_pos += 1
                    text_pos += 1
                elif op == 'I':
                    text_pos += 1
                elif op == 'D':
                    pattern_pos += 1
                else:
                    print(f"Invalid op {op}")
                    return False
    except IndexError:
        # Reading outside the pattern or text
        return False

    if (pattern_pos != len(pattern)) or (text_pos != len(text)):
        return False

    return True

def check_results_file(results_file, lines, config):
    """Validate one results file without writing to shared terminal state."""
    result = FileCheckResult(results_file)
    if lines is None:
        try:
            with open(results_file, 'r') as results_handle:
                lines = results_handle.readlines()
        except (OSError, UnicodeError) as error:
            result.errors.append(f"Error opening file {results_file}: {error}")
            result.return_code = 1
            return result

    gt_scores = config.ground_truth_scores
    sequence_pairs = config.sequence_pairs
    if gt_scores is not None and len(gt_scores) != len(lines):
        result.errors.append(
            f"Ground truth has {len(gt_scores)} rows but {results_file} has "
            f"{len(lines)} rows"
        )
        result.return_code = 1
    if sequence_pairs is not None and len(sequence_pairs) != len(lines):
        result.errors.append(
            f"Sequences file has {len(sequence_pairs)} pairs but {results_file} "
            f"has {len(lines)} rows"
        )
        result.return_code = 1

    for line_num, line in enumerate(lines):
        elements = line.rstrip().split()
        if len(elements) not in (2, 4):
            result.errors.append(f"Invalid result row at line {line_num + 1}")
            result.incorrect += 1
            result.return_code = 1
            continue

        try:
            score = abs(int(elements[0]))
            cigar = elements[1]
        except ValueError:
            result.errors.append(f"Invalid score at line {line_num + 1}")
            result.incorrect += 1
            result.return_code = 1
            continue

        if config.collect_plot_scores:
            result.plot_scores.append(score)

        try:
            ops, cigar_reps = parse_cigar(cigar)
        except ValueError as error:
            result.errors.append(
                f"Invalid CIGAR at line {line_num + 1}: {error}"
            )
            result.incorrect += 1
            result.return_code = 1
            continue

        is_traceback_correct = True
        if len(elements) == 4:
            pattern, text = elements[2:4]
            is_traceback_correct = check_cigar_sequences(
                score, ops, cigar_reps, pattern, text, config.with_mismatches
            )
        elif sequence_pairs is not None:
            if line_num >= len(sequence_pairs):
                is_traceback_correct = False
            else:
                pattern, text = sequence_pairs[line_num]
                is_traceback_correct = check_cigar_sequences(
                    score, ops, cigar_reps, pattern, text,
                    config.with_mismatches
                )

        if config.distance_function == 'edit':
            is_score_correct, cigar_score = check_score_edit(
                score, ops, cigar_reps
            )
        elif config.distance_function == 'affine':
            is_score_correct, cigar_score = check_score_affine(
                score, ops, cigar_reps, *config.penalties[:4]
            )
        else:
            is_score_correct, cigar_score = check_score_affine2p(
                score, ops, cigar_reps, *config.penalties[:6]
            )

        is_ground_truth_correct = (
            gt_scores is None
            or (line_num < len(gt_scores) and gt_scores[line_num] == score)
        )
        is_correct = (
            is_score_correct
            and is_traceback_correct
            and is_ground_truth_correct
        )

        if is_correct:
            result.correct += 1
        else:
            gt_score = (
                gt_scores[line_num]
                if gt_scores is not None and line_num < len(gt_scores)
                else None
            )
            result.incorrect_rows.append(
                (line_num + 1, score, cigar, cigar_score, gt_score)
            )
            result.incorrect += 1

        result.score_total += score
        result.score_count += 1
        result.max_score = max(result.max_score, score)

    if result.correct + result.incorrect == 0:
        result.errors.append(f"No valid CIGARs found in {results_file}")
    if result.incorrect > 0 or result.score_count == 0:
        result.return_code = 1
    return result


def initialize_worker(config):
    global _WORKER_CONFIG
    _WORKER_CONFIG = config


def check_results_file_in_worker(results_file, lines):
    return check_results_file(results_file, lines, _WORKER_CONFIG)


def run_file_checks(files, stdin_lines, config, jobs, quiet):
    """Run file checks, using one parent-owned progress bar."""
    inputs = [
        (filename, stdin_lines if filename == '-' else None, config)
        for filename in files
    ]
    completed = [None] * len(inputs)
    progress = tqdm(
        total=len(inputs), unit='files', ncols=120, disable=quiet, leave=False,
        desc='Checking files'
    )
    if jobs == 1 or len(inputs) == 1:
        for index, check_input in enumerate(inputs):
            result = check_results_file(*check_input)
            completed[index] = result
            progress.set_postfix_str(os.path.basename(result.filename))
            progress.update()
    else:
        start_method = (
            'fork' if 'fork' in multiprocessing.get_all_start_methods()
            else 'spawn'
        )
        context = multiprocessing.get_context(start_method)
        with ProcessPoolExecutor(
            max_workers=min(jobs, len(inputs)), mp_context=context,
            initializer=initialize_worker, initargs=(config,)
        ) as executor:
            futures = {
                executor.submit(check_results_file_in_worker, *check_input[:2]):
                    (index, check_input[0])
                for index, check_input in enumerate(inputs)
            }
            for future in as_completed(futures):
                index, filename = futures[future]
                try:
                    result = future.result()
                except Exception as error:
                    result = FileCheckResult(
                        filename, return_code=1,
                        errors=[f"Error checking file {filename}: {error}"]
                    )
                completed[index] = result
                progress.set_postfix_str(os.path.basename(filename))
                progress.update()
    progress.close()
    return completed


def checkalign():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='*', help='Files with the results to check (- for stdin)')
    parser.add_argument('-g', '--penalties', default='0,1,0,1,0,0', help='Penalties in a,x,o,e,o1,e1 format (match, mismatch, gap-open, gap-extend, gap-open1, gap-extend1). Default is 0,1,0,1,0,0 (equivalent to edit distance)')
    parser.add_argument('-d', '--distance-function', choices=('edit', 'affine', 'affine2p'), default='edit', help='Distance function. \'edit\', \'affine\' or \'affine2p\'. Default is \'edit\'')
    parser.add_argument('-j', '--jobs', type=positive_integer, default=1, help='Number of result files to check in parallel. Default is 1')
    parser.add_argument('-q', '--quiet', required=False, action='store_true', help='Don\'t print any output on the stdout')
    parser.add_argument('-v', '--verbose', required=False, action='store_true', help='Print additonal information about incorrect CIGARs')
    parser.add_argument('-s', '--sequences', required=False, help='File with the input sequences')
    parser.add_argument('-t', '--ground-truth', required=False, help='File with the ground truth')
    parser.add_argument('-x', '--no-mismatches', required=False, action='store_true', help='Check CIGARs without mismatches (e.g. the ones produced by KSW2)')
    parser.add_argument('-p', '--plot', required=False, action='store_true', help='Create a plot with cumulative score')
    args = parser.parse_args()

    try:
        penalty_values = list(map(int, args.penalties.split(',')))
    except ValueError:
        parser.error("Penalties must be comma-separated integers")
    if args.distance_function == 'affine' and len(penalty_values) < 4:
        parser.error("Affine distance requires at least four penalties")
    if args.distance_function == 'affine2p' and len(penalty_values) < 6:
        parser.error("Affine2p distance requires at least six penalties")
    if args.distance_function in ('affine', 'affine2p'):
        required = 4 if args.distance_function == 'affine' else 6
        penalty_values = penalty_values[:required]
        sign = -1 if penalty_values[0] >= 0 else 1
        penalty_values[1:] = [sign * abs(value)
                              for value in penalty_values[1:]]

    if not args.files:
        parser.print_help()
        sys.exit(1)
    if args.files.count('-') > 1:
        parser.error("Standard input (-) can only be specified once")

    sequence_pairs = None
    if args.sequences:
        try:
            with open(args.sequences, 'r') as seq_handle:
                sequence_pairs = tuple(parse_sequence_pairs(seq_handle.readlines()))
        except (OSError, UnicodeError, ValueError) as error:
            error_console.print(f"Error reading sequences file: {error}")
            sys.exit(1)

    gt_scores = None
    if args.ground_truth:
        gt_scores = []
        try:
            with open(args.ground_truth, 'r') as gt_handle:
                for line_num, line in enumerate(gt_handle, start=1):
                    try:
                        gt_scores.append(abs(int(line.split()[0])))
                    except (IndexError, ValueError):
                        error_console.print(
                            f"Invalid ground-truth score at line {line_num}"
                        )
                        sys.exit(1)
        except (OSError, UnicodeError) as error:
            error_console.print(
                f"Error opening ground-truth file {args.ground_truth}: {error}"
            )
            sys.exit(1)
        gt_scores = tuple(gt_scores)

    config = CheckConfig(
        distance_function=args.distance_function,
        penalties=tuple(penalty_values),
        with_mismatches=not args.no_mismatches,
        ground_truth_scores=gt_scores,
        sequence_pairs=sequence_pairs,
        collect_plot_scores=args.plot,
    )
    stdin_lines = sys.stdin.readlines() if '-' in args.files else None
    file_results = run_file_checks(
        args.files, stdin_lines, config, args.jobs, args.quiet
    )

    results_table = Table(title="Results")
    results_table.add_column("File")
    results_table.add_column("Correct", style="green")
    results_table.add_column("Incorrect", style="red")
    results_table.add_column("Accuracy", style="bold")
    plot_data = {'ground_truth': list(gt_scores or ())} if args.plot else None
    return_code = 0

    for result in file_results:
        return_code = max(return_code, result.return_code)
        for error in result.errors:
            error_console.print(error)
        total = result.correct + result.incorrect
        if total:
            results_table.add_row(
                f'{"✅" if result.incorrect == 0 else "❌"} {result.filename}',
                str(result.correct), str(result.incorrect),
                f"{result.correct * 100 / total:.2f}%"
            )
        if result.incorrect and args.verbose and not args.quiet:
            incorrect_table = generate_incorrect_cigars_table(
                result.filename, gt_scores is not None
            )
            for incorrect_row in result.incorrect_rows:
                update_incorrect_cigars_table(incorrect_table, *incorrect_row)
            console.print(incorrect_table)
        if not args.quiet and result.score_count:
            print(
                f'Average score for {result.filename}: '
                f'{result.score_total / result.score_count:.2f}'
            )
            print(f'Max score for {result.filename}: {result.max_score}')
        if args.plot:
            plot_data[result.filename] = result.plot_scores

    if args.plot:
        plot_cummulative_scores(plot_data)
    if not args.quiet:
        console.print(results_table)
    sys.exit(return_code)

if __name__ == '__main__':
    checkalign()
