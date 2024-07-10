import argparse
import numpy as np
import matplotlib.pyplot as plt

def show_cigar(cigar):
    tlen = 0
    plen = 0

    # Get the repetitions of each operation
    tmp_cigar = cigar.replace('M', ' ').replace('I', ' ').replace('D', ' ').replace('X', ' ')
    try:
        cigar_reps = list(map(int, tmp_cigar.split()))
    except ValueError:
        print(f"Invalid CIGAR")
        return

    # Get the operations
    ops = []
    for e in cigar:
        if e in ['M', 'X', 'I', 'D']:
            ops.append(e) 

    score = 0
    # Get tlen, plen, and score
    for i in range(len(cigar_reps)):
        curr_rep = cigar_reps[i]
        curr_op = ops[i]
        for _ in range(curr_rep):
            # Set plen (vertical) and tlen (horizontal)
            if curr_op == 'M':
                plen += 1
                tlen += 1
            elif curr_op == 'X':
                plen += 1
                tlen += 1
                score += 1
            elif curr_op == 'I':
                plen += 1
                score += 1
            elif curr_op == 'D':
                tlen += 1
                score += 1

    v, h = 0, 0
    coords = []
    vals = []

    partial_score = 0
    bps = []
    bp_cutoff = score // 10

    for i in range(len(cigar_reps)):
        curr_rep = cigar_reps[i]
        curr_op = ops[i]
        for _ in range(curr_rep):
            try:
                if curr_op == 'X':
                    coords.append((h, v))
                    vals.append(2)
                    v += 1
                    h += 1
                    partial_score += 1
                elif curr_op == 'M':
                    coords.append((h, v))
                    vals.append(1)
                    v += 1
                    h += 1
                elif curr_op == 'I':
                    coords.append((h, v))
                    vals.append(3)
                    v += 1
                    partial_score += 1
                elif curr_op == 'D':
                    coords.append((h, v))
                    vals.append(4)
                    h += 1
                    partial_score += 1

                if partial_score % bp_cutoff == 0:
                    bps.append((h, v))
            except IndexError:
                print(f"IndexError: v={v} h={h}. plen={plen} tlen={tlen}")
                return


    print('Matrix constructed')

    aspect_ratio = tlen/plen
    fig, ax = plt.subplots(figsize=(30, 30/aspect_ratio))

    #cmap = plt.cm.viridis
    #cmap.set_under(color='white')
    #cmap.set_over(color='red')
    #ax.imshow(m, vmin=0.1, vmax=4, cmap=cmap)

    xs = [coord[0] for coord in coords]
    ys = [coord[1] for coord in coords]

    ax.set_xlim(0, tlen)
    ax.set_ylim(0, plen)
    ax.xaxis.tick_top()

    plt.gca().invert_yaxis()

    # colors, vals==1: red, vals==2: blue, vals==3: green, vals==4: yellow
    c = [_, 'red', 'blue', 'green', 'yellow']
    colors = [c[val] for val in vals]

    ax.scatter(xs, ys, s=0.001, c=colors)

    for i in range(0, tlen, int(tlen*0.1)):
        ax.axhline(i, color='lightgray', lw=1, ls='--')
    for i in range(0, plen, int(plen*0.1)):
        ax.axvline(i, color='lightgray', lw=1, ls='--')

    # Legend: red=M, blue=X, green=I, yellow=D
    ax.plot([0], [0], 'ro', label='M')
    ax.plot([0], [0], 'bo', label='X')
    ax.plot([0], [0], 'go', label='I')
    ax.plot([0], [0], 'yo', label='D')
    ax.legend()

    # Print score in the figure
    ax.text(0, plen, f"Score: {score}", fontsize=24, color='black')

    # Print read line from origin (0, 0) to end (tlen, plen)
    ax.plot([0, tlen], [0, plen], color='red', lw=1, ls='--')

    # Plot breakpoints
    for bp in bps:
        ax.plot(bp[0], bp[1], 'ro', markersize=5)

    # tight layout
    plt.tight_layout()

    fig.savefig('cigar.png')
    plt.show()

def catcigar():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs='?', help='Input file in .out format')
    parser.add_argument('-i', '--input', required=True, help='Input file in .out format')
    parser.add_argument('-n', '--number', type=int, help='Number of CIGARS to read')
    parser.add_argument('-s', '--skip', type=int, help='Number of CIGARS to skip')

    args = parser.parse_args()

    with open(args.input, 'r') as f:
        count = 0
        for idx, line in enumerate(f):
            if args.skip and idx < args.skip:
                continue
            if args.number and count == args.number:
                break
            score, cigar = line.strip().split()[:2]
            print(f"Score: {score} CIGAR: {cigar}")
            show_cigar(cigar)
            count += 1

if __name__ == '__main__':
    catcigar()