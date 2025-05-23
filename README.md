# PAtools

Collection of tools and scripts to work with Pairwise Alignment data, specially with `.seq` and `.out` files. Still work in progress.

## Installation

Clone the repository:
Optionally, create a virtual environment to avoid conflicts with other Python packages:
```bash
python3 -m venv venv
source venv/bin/activate
```

Install the tools:
```bash
pip3 install git+https://github.com/quim0/patools.git
```

## Tools

### `checkalign`

Check for the correctness of one or multiple `.out` files, or from stdin.

```
usage: checkalign.py [-h] [-g PENALTIES] [-d DISTANCE_FUNCTION] [-q] [-v]
                     [-s SEQUENCES] [-t GROUND_TRUTH] [-p]
                     [files ...]

positional arguments:
  files                 Files with the results to check (- for stdin)

options:
  -h, --help            show this help message and exit
  -g PENALTIES, --penalties PENALTIES
                        Penalties in a,x,o,e,o1,e1 format (match, mismatch,
                        gap-open, gap-extend, gap-open1, gap-extend1). Default
                        is 0,1,0,1,0,0 (equivalent to edit distance)
  -d DISTANCE_FUNCTION, --distance-function DISTANCE_FUNCTION
                        Distance function. 'edit', 'affine' or 'affine2p'.
                        Default is 'edit'
  -q, --quiet           Don't print any output on the stdout
  -v, --verbose         Print additonal information about incorrect CIGARs
  -s SEQUENCES, --sequences SEQUENCES
                        File with the input sequences
  -t GROUND_TRUTH, --ground-truth GROUND_TRUTH
                        File with the ground truth
  -p, --plot            Create a plot with cumulative score
```

`.out` file format, where `SEQUENCE1` and `SEQUENCE2` are optional, and spaces can be replaced by tabs:
```
SCORE CIGAR SEQUENCE1 SEQUENCE2
SCORE CIGAR SEQUENCE1 SEQUENCE2
SCORE CIGAR SEQUENCE1 SEQUENCE2
```

### `catcigar`

Takes one or multiple CIGARS from an `.out` file and generates an illustration of the alignment.

```
usage: catcigar [-h] -i INPUT [-n NUMBER] [-s SKIP] [file]

positional arguments:
  file                  Input file in .out format

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file in .out format
  -n NUMBER, --number NUMBER
                        Number of CIGARS to read
  -s SKIP, --skip SKIP  Number of CIGARS to skip
```

## Examples

Check if a file `cigars.out` is correct, and compare them with a ground truth `correct_cigars.out`:
```bash
checkalign -t correct_cigars.out cigars.out
```