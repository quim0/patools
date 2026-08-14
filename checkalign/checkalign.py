#!/usr/bin/env python3
import argparse
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

def checkalign():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='*', help='Files with the results to check (- for stdin)')
    parser.add_argument('-g', '--penalties', default='0,1,0,1,0,0', help='Penalties in a,x,o,e,o1,e1 format (match, mismatch, gap-open, gap-extend, gap-open1, gap-extend1). Default is 0,1,0,1,0,0 (equivalent to edit distance)')
    parser.add_argument('-d', '--distance-function', choices=('edit', 'affine', 'affine2p'), default='edit', help='Distance function. \'edit\', \'affine\' or \'affine2p\'. Default is \'edit\'')
    parser.add_argument('-q', '--quiet', required=False, action='store_true', help='Don\'t print any output on the stdout')
    parser.add_argument('-v', '--verbose', required=False, action='store_true', help='Print additonal information about incorrect CIGARs')
    parser.add_argument('-s', '--sequences', required=False, help='File with the input sequences')
    parser.add_argument('-t', '--ground-truth', required=False, help='File with the ground truth')
    parser.add_argument('-x', '--no-mismatches', required=False, action='store_true', help='Check CIGARs without mismatches (e.g. the ones produced by KSW2)')
    parser.add_argument('-p', '--plot', required=False, action='store_true', help='Create a plot with cumulative score')

    args = parser.parse_args()

    penalties = args.penalties.split(',')
    try:
        penalty_values = list(map(int, penalties))
    except ValueError:
        parser.error("Penalties must be comma-separated integers")

    if args.distance_function == 'affine' and len(penalties) < 4:
        parser.error("Affine distance requires at least four penalties")
    if args.distance_function == 'affine2p' and len(penalties) < 6:
        parser.error("Affine2p distance requires at least six penalties")

    if args.distance_function == 'affine':
        M,X,O,E = penalty_values[:4]
        if M >= 0: 
            X,O,E = -abs(X),-abs(O),-abs(E)
        else:
            X,O,E = abs(X),abs(O),abs(E)
    elif args.distance_function == 'affine2p':
        M,X,O,E,O1,E1 = penalty_values[:6]
        if M >= 0: 
            X,O,E,O1,E1 = -abs(X),-abs(O),-abs(E),-abs(O1),-abs(E1)
        else:
            X,O,E,O1,E1 = abs(X),abs(O),abs(E),abs(O1),abs(E1)

    if args.files == []:
        parser.print_help()
        sys.exit(1)

    sequence_pairs = None
    if args.sequences:
        try:
            with open(args.sequences, 'r') as seq_f:
                sequence_lines = seq_f.readlines()
        except (OSError, UnicodeError) as e:
            error_console.print(f"Error opening sequences file: {e}")
            sys.exit(1)

        try:
            sequence_pairs = parse_sequence_pairs(sequence_lines)
        except ValueError as e:
            error_console.print(str(e))
            sys.exit(1)

    plot_data = {}
    if args.plot:
        # For each file, create dict entry to store all the scores and another for
        # the ground truth.
        plot_data = {f: [] for f in args.files}
        plot_data['ground_truth'] = []

    with_mismatches = True
    if args.no_mismatches:
        with_mismatches = False

    with_ground_truth = False
    gt_scores = []
    if args.ground_truth:
        curr_gt = args.ground_truth
        try:
            with open(curr_gt, 'r') as fgt:
                for gt_line_num, line in enumerate(fgt, start=1):
                    elements = line.split()
                    try:
                        gt_scores.append(abs(int(elements[0])))
                    except (IndexError, ValueError):
                        error_console.print(
                            f"Invalid ground-truth score at line {gt_line_num}"
                        )
                        sys.exit(1)
        except (OSError, UnicodeError) as e:
            error_console.print(f'Error opening ground-truth file {curr_gt}: {e}')
            sys.exit(1)

        with_ground_truth = True
        if args.plot:
            plot_data['ground_truth'] = gt_scores

    retval = 0
    results = Table(title="Results")
    results.add_column("File")
    results.add_column("Correct", style="green")
    results.add_column("Incorrect", style="red")
    results.add_column("Accuracy", style="bold")

    for f in args.files:
        results_file = f
        if results_file == '-':
            lines = sys.stdin.readlines()
        else:
            try:
                with open(results_file, 'r') as f:
                    lines = f.readlines()
            except (FileNotFoundError, IsADirectoryError) as e:
                error_console.print(f'Error opening file {results_file}... Skipping.')
                sys.exit(1)

        score_total = 0
        score_count = 0
        max_score = 0
        correct = 0
        incorrect = 0
        incorrect_cigars_table = generate_incorrect_cigars_table(results_file, with_ground_truth)
        pbar = tqdm(total=len(lines), unit='CIGARs', ncols=120,
                    disable=args.quiet, leave=False,
                    bar_format='{l_bar}{bar}{r_bar}' + f' {os.path.basename(results_file)}')
        if with_ground_truth and len(gt_scores) != len(lines):
            error_console.print(
                f"Ground truth has {len(gt_scores)} rows but {results_file} has "
                f"{len(lines)} rows"
            )
            retval = 1
        if sequence_pairs is not None and len(sequence_pairs) != len(lines):
            error_console.print(
                f"Sequences file has {len(sequence_pairs)} pairs but {results_file} "
                f"has {len(lines)} rows"
            )
            retval = 1

        for line_num, line in enumerate(lines):
            pbar.set_description(f'(correct={correct}, incorrect={incorrect})')
            line = line.rstrip()
            elements = line.split()
            if len(elements) not in (2, 4):
                error_console.print(f"Invalid result row at line {line_num + 1}")
                incorrect += 1
                retval = 1
                pbar.update(1)
                continue

            try:
                score = abs(int(elements[0]))
                cigar = elements[1]
            except ValueError:
                error_console.print(f"Invalid score at line {line_num + 1}")
                incorrect += 1
                retval = 1
                pbar.update(1)
                continue

            if args.plot:
                plot_data[results_file].append(score)

            try:
                ops, cigar_reps = parse_cigar(cigar)
            except ValueError as e:
                error_console.print(f"Invalid CIGAR at line {line_num + 1}: {e}")
                incorrect += 1
                retval = 1
                pbar.update(1)
                continue

            is_traceback_correct = True
            if len(elements) == 4:
                pattern = elements[2]
                text = elements[3]
                is_traceback_correct = check_cigar_sequences(
                    score, ops, cigar_reps, pattern, text, with_mismatches
                )
            elif sequence_pairs is not None:
                if line_num >= len(sequence_pairs):
                    is_traceback_correct = False
                else:
                    pattern, text = sequence_pairs[line_num]
                    is_traceback_correct = check_cigar_sequences(
                        score, ops, cigar_reps, pattern, text, with_mismatches
                    )

            if args.distance_function == 'edit':
                is_score_correct, cigar_score = check_score_edit(score, ops, cigar_reps)
            elif args.distance_function == 'affine':
                is_score_correct, cigar_score = check_score_affine(score, ops, cigar_reps, M, X, O, E)
            else:
                is_score_correct, cigar_score = check_score_affine2p(
                    score, ops, cigar_reps, M, X, O, E, O1, E1
                )

            is_ground_truth_correct = (
                not with_ground_truth
                or (line_num < len(gt_scores) and gt_scores[line_num] == score)
            )
            is_correct = (
                is_score_correct
                and is_traceback_correct
                and is_ground_truth_correct
            )

            if not is_correct or not is_traceback_correct:
                gt_score = (
                    gt_scores[line_num]
                    if with_ground_truth and line_num < len(gt_scores)
                    else None
                )
                update_incorrect_cigars_table(
                    incorrect_cigars_table, line_num + 1, score, cigar,
                    cigar_score, gt_score
                )
                #if not args.quiet:
                #    get_incorrect_cigar_table(line_num, score, cigar, cigar_score, gt_scores[line_num] if with_ground_truth else None)
                incorrect += 1
            else:
                correct += 1

            score_total += score
            score_count += 1
            max_score = max(max_score, score)

            pbar.update(1)
        pbar.close()

        if (correct+incorrect) == 0:
            error_console.print(f"[red]No valid CIGARs found in {results_file}[/red]")
        else:
            results.add_row(f'{"✅" if incorrect==0 else "❌"} {results_file}', str(correct), str(incorrect), f"{correct*100/(correct+incorrect):.2f}%")

        if incorrect > 0 and args.verbose and not args.quiet:
            console.print(incorrect_cigars_table)

        if not args.quiet and score_count:
            print(f'Average score for {results_file}: {score_total/score_count:.2f}')
            print(f'Max score for {results_file}: {max_score}')

        if incorrect > 0 or score_count == 0:
            retval = 1

    if args.plot:
        plot_cummulative_scores(plot_data)

    if not args.quiet:
        console.print(results)

    sys.exit(retval)

if __name__ == '__main__':
    checkalign()
