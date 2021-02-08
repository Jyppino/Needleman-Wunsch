# Needleman-Wunsch

Basic Python script performing the [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) for sequence alignment.
Returns max similarity matrix and either upmost or downmost alignment of two sequences.

## Usage

```
usage: python wunsch.py [-h] [--gap [GAP]] [--match [MATCH]] [--mismatch [MISMATCH]] [--down] seq_1 seq_2

Needleman-Wunsch Script

positional arguments:
  seq_1                 String of first sequence
  seq_2                 String of second sequence

optional arguments:
  -h, --help            show this help message and exit
  --gap [GAP]           Cost of gaps, default=-2
  --match [MATCH]       Reward of match, default=1
  --mismatch [MISMATCH]
                        Cost of mismatch, default=-1
  --down                Use downmost alignment instead of upmost, default=False
```

## Example

`python wunsch.py 'ACCTGGGCCACGT' 'ACCAGGCTACGA'`

Max similarity matrix:

```
[[  0  -2  -4  -6  -8 -10 -12 -14 -16 -18 -20 -22 -24 -26]
 [ -2   1  -1  -3  -5  -7  -9 -11 -13 -15 -17 -19 -21 -23]
 [ -4  -1   2   0  -2  -4  -6  -8 -10 -12 -14 -16 -18 -20]
 [ -6  -3   0   3   1  -1  -3  -5  -7  -9 -11 -13 -15 -17]
 [ -8  -5  -2   1   2   0  -2  -4  -6  -8  -8 -10 -12 -14]
 [-10  -7  -4  -1   0   3   1  -1  -3  -5  -7  -9  -9 -11]
 [-12  -9  -6  -3  -2   1   4   2   0  -2  -4  -6  -8 -10]
 [-14 -11  -8  -5  -4  -1   2   3   3   1  -1  -3  -5  -7]
 [-16 -13 -10  -7  -4  -3   0   1   2   2   0  -2  -4  -4]
 [-18 -15 -12  -9  -6  -5  -2  -1   0   1   3   1  -1  -3]
 [-20 -17 -14 -11  -8  -7  -4  -3   0   1   1   4   2   0]
 [-22 -19 -16 -13 -10  -7  -6  -3  -2  -1   0   2   5   3]
 [-24 -21 -18 -15 -12  -9  -8  -5  -4  -3   0   0   3   4]]
```

Upmost alignment:

```
ACC-AGGCTACGA
ACCTGGGCCACGT
```

Downmost alignment:

```
ACCAGG-CTACGA
ACCTGGGCCACGT
```
