#!/usr/bin/env python3
import argparse
import sys, os
from tqdm import tqdm
from rich.console import Console
from rich.table import Table

PROGRESS_CHECK = 5000

console = Console()
error_console = Console(stderr=True)

def print_report(correct, incorrect, results_file):
    console.print(f"([bold]{results_file}[/bold]) Correct=[green]{correct}[/green], Incorrect=[red]{incorrect}[/red], Accuracy=[bold]{correct*100/(correct+incorrect):.2f}%[/bold]")

def update_incorrect_cigars_table(table, line_num, score, cigar, cigar_score, gt_score=None):
    args = [str(line_num), str(score), cigar, str(cigar_score)]
    if gt_score:
        args.append(str(gt_score))
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

def check_score_affine(score, cigar_ops, cigar_reps, X, O, E):
    score_calc = 0
    for idx, op in enumerate(cigar_ops):
        reps = cigar_reps[idx]
        if op == 'M':
            continue
        elif op == 'X':
            score_calc += X * reps
        elif op in ['I', 'D']:
            score_calc += O + E * reps

    return (score == score_calc, score_calc)

def check_score_affine2p(score, cigar_ops, cigar_reps, X, O1, E1, O2, E2):
    score_calc = 0
    for idx, op in enumerate(cigar_ops):
        reps = cigar_reps[idx]
        if op == 'M':
            continue
        elif op == 'X':
            score_calc += X * reps
        elif op in ['I', 'D']:
            score_calc += min(
                O1 + E1 * reps,
                O2 + E2 * reps
            )

    return (score == score_calc, score_calc)

def check_cigar_sequences(score, cigar_ops, cigar_reps, pattern, text):
    text_pos = 0
    pattern_pos = 0

    try:
        for idx, op in enumerate(cigar_ops):
            reps = cigar_reps[idx]
            for _ in range(reps):
                if op == 'M':
                    if pattern[pattern_pos] != text[text_pos]:
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
    parser.add_argument('-g', '--penalties', default='1,0,1,0,0', help='Penalties in x,o,e,o1,e1 format (mismatch, gap-open, gap-extend, gap-open1, gap-extend1). Default is 1,0,1,0,0 (equivalent to edit distance)')
    parser.add_argument('-d', '--distance-function', default='edit', help='Distance function. \'edit\', \'affine\' or \'affine2p\'. Default is \'edit\'')
    parser.add_argument('-q', '--quiet', required=False, action='store_true', help='Don\'t print any output on the stdout')
    parser.add_argument('-v', '--verbose', required=False, action='store_true', help='Print additonal information about incorrect CIGARs')
    parser.add_argument('-s', '--sequences', required=False, help='File with the input sequences')
    parser.add_argument('-t', '--ground-truth', required=False, help='File with the ground truth')
    parser.add_argument('-p', '--plot', required=False, action='store_true', help='Create a plot with cumulative score')

    args = parser.parse_args()

    penalties = args.penalties.split(',')
    if args.distance_function == 'affine' and len(penalties) < 3:
        error_console.print("Invalid number of penalties")
        quit(1)
    if args.distance_function == 'affine2p' and len(penalties) < 5:
        error_console.print("Invalid number of penalties for affine2p")
        quit(1)

    if args.distance_function == 'affine':
        X,O,E = map(int, penalties[:3])
    elif args.distance_function == 'affine2p':
        X,O,E,O1,E1 = map(int, penalties[:5])

    if args.files == []:
        parser.print_help()
        quit(1)

    # Open sequences file on the fly if needed
    seq_f = None
    if args.sequences:
        try:
            seq_f = open(args.sequences)
        except Exception as e:
            print(f"Error opening sequences file: {e}")
            seq_f = None

    plot_data = {}
    if args.plot:
        # For each file, create dict entry to store all the scores and another for
        # the ground truth.
        plot_data = {f: [] for f in args.files}
        for f in args.files:
            plot_data[f'ground_truth'] = []

    with_ground_truth = False
    if args.ground_truth:
        curr_gt = args.ground_truth
        if not curr_gt:
            error_console.print(f"No ground truth for {f}.")

        else:
            with_ground_truth = True
            gt_scores = []
            try:
                with open(curr_gt, 'r') as fgt:
                    for l in fgt:
                        l = l.rstrip()
                        elements = l.split()
                        try:
                            score = abs(int(elements[0]))
                        except ValueError:
                            print(f"Invalid score at line {line_num}")
                            exit(1)
                        gt_scores.append(score)

                    # If plotting in enabled, store the ground truth
                    if args.plot:
                        plot_data['ground_truth'] = gt_scores

            except (FileNotFoundError, IsADirectoryError) as e:
                error_console.print(f'Error opening file {curr_gt}... Skipping.')
                with_ground_truth = False

    reval = 0
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
                quit(1)

        avg_score = 0
        max_score = 0
        correct = 0
        incorrect = 0
        incorrect_cigars_table = generate_incorrect_cigars_table(results_file, with_ground_truth)
        pbar = tqdm(total=len(lines), unit='CIGARs', ncols=120,
                    disable=args.quiet, leave=False,
                    bar_format='{l_bar}{bar}{r_bar}' + f' {os.path.basename(results_file)}')
        for line_num, line in enumerate(lines):
            pbar.set_description(f'(correct={correct}, incorrect={incorrect})')
            line = line.rstrip()
            elements = line.split()
            if len(elements) < 2 or len(elements) > 4:
                error_console.print(f"Invalid score or CIGAR at line {line_num}")
                pbar.update(1)
                continue

            try:
                score = abs(int(elements[0]))
                cigar = elements[1]
            except ValueError:
                error_console.print(f"Invalid score or CIGAR at line {line_num}")
                pbar.update(1)
                continue

            if args.plot:
                plot_data[results_file].append(score)

            cigar_tmp = cigar.replace('M', ' ')
            cigar_tmp = cigar_tmp.replace('X', ' ')
            cigar_tmp = cigar_tmp.replace('I', ' ')
            cigar_tmp = cigar_tmp.replace('D', ' ')
            try:
                cigar_reps = list(map(int, cigar_tmp.split()))
            except ValueError:
                error_console.print(f"Invalid op at CIGAR {line_num}")
                pbar.update(1)
                continue

            ops = []
            for e in cigar:
                if e in ['M', 'X', 'I', 'D']:
                    ops.append(e)

            is_correct = True

            if with_ground_truth:
                if gt_scores[line_num] != score:
                    is_correct = False

            if len(elements) == 4:
                pattern = elements[2]
                text = elements[3]
                ok = check_cigar_sequences(score, ops, cigar_reps, pattern, text)
                if not ok:
                    is_correct = False
                    #if not args.quiet:
                    #    error_console.print(f"CIGAR {line_num} do not fit the pattern and text.")
            elif seq_f:
                # Read next two lines from seq_f
                pat_line = seq_f.readline()
                txt_line = seq_f.readline()
                pattern = text = ''
                if not pat_line or not txt_line:
                    print(f"Sequences file ended prematurely at line {line_num}")
                    seq_f.close()
                    seq_f = None
                    pattern = text = ''
                    break
                else:
                    pattern = pat_line.strip()[1:]
                    text    = txt_line.strip()[1:]
                ok = check_cigar_sequences(score, ops, cigar_reps, pattern, text)
                if not ok:
                    is_correct = False
                    #if not args.quiet:
                    #    error_console.print(f"CIGAR {line_num} do not fit the pattern and text.")

            # Calculate score
            if args.distance_function == 'edit':
                is_correct, cigar_score = check_score_edit(score, ops, cigar_reps)
            elif args.distance_function == 'affine':
                is_correct, cigar_score = check_score_affine(score, ops, cigar_reps, X, O, E)
            elif args.distance_function == 'affine2p':
                is_correct, cigar_score = check_score_affine2p(score, ops, cigar_reps, X, O, E, O1, E1)
            else:
                error_console.print(f"Invalid distance function {args.distance_function}")
                quit

            if not is_correct:
                update_incorrect_cigars_table(incorrect_cigars_table, line_num, score, cigar, cigar_score, gt_scores[line_num] if with_ground_truth else None)
                #if not args.quiet:
                #    get_incorrect_cigar_table(line_num, score, cigar, cigar_score, gt_scores[line_num] if with_ground_truth else None)
                incorrect += 1
            else:
                correct += 1

            avg_score += score
            max_score = max(max_score, score)

            pbar.update(1)

        pbar.close()

        if (correct+incorrect) == 0:
            error_console.print(f"[red]No valid CIGARs found in {results_file}[/red]")
        else:
            results.add_row(f'{"✅" if incorrect==0 else "❌"} {results_file}', str(correct), str(incorrect), f"{correct*100/(correct+incorrect):.2f}%")

        if incorrect > 0 and args.verbose and not args.quiet:
            console.print(incorrect_cigars_table)

        print(f'Average score for {results_file}: {avg_score/(correct+incorrect):.2f}')
        print(f'Max score for {results_file}: {max_score}')

        if args.plot:
            plot_cummulative_scores(plot_data)

        if incorrect > 0:
            retval = 1

    if not args.quiet:
        console.print(results)

    quit(reval)

if __name__ == '__main__':
    checkalign()