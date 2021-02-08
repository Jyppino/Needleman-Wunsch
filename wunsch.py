# Needleman-Wunsch Script
# Jaap de Boer
# UCD 20203421
# ./wunsch -h for usage information

import numpy as np
import argparse

# Main method
def main(argv):
    global SEQ1, SEQ2, GAP, MATCH, MISMATCH

    # Setup global variables (constants)
    SEQ1 = argv.seq_1
    SEQ2 = argv.seq_2
    GAP = argv.gap
    MATCH = argv.match
    MISMATCH = argv.mismatch

    a = wunsch()
    
    align1, align2 = [], []
    if not argv.down:
        (align1, align2) = alignUpMost(len(SEQ1), len(SEQ2), a)
    else:
        (align1, align2) = alignDownMost(len(SEQ1), len(SEQ2), a)

    print(a)
    print(''.join(align1))
    print(''.join(align2))

# Create max similarity matrix (Needleman-Wunsch algorithm)
# Returns similarity matrix
def wunsch():
    x = len(SEQ1)
    y = len(SEQ2)
    a = np.zeros((y+1, x+1), dtype=int)

    for i in range(x+1):
        a[0][i] = i * GAP
    for j in range(y+1):
        a[j][0] = j * GAP
    
    for i in range(1, y+1):
        for j in range(1, x+1):
            a[i][j] = max(
                a[i-1][j] + GAP,
                a[i][j-1] + GAP,
                a[i-1][j-1] + calcPenalty(j-1, i-1)
            )

    return a

# Calculate penalty
# x, y: positions in SEQ1 and SEQ2
# Returns integer (score)
def calcPenalty(x, y):
    return MATCH if SEQ1[x] == SEQ2[y] else MISMATCH

# Backtrack steps for best global upmost alignments
# x, y: current coordinates (recursion)
# a: similarity matrix
# returns tuple of lists (alignments)
def alignUpMost(x, y, a, allY=[], allX=[]):
    if x == 0 and y == 0:
        return (allY, allX)
    elif y > 0 and a[y][x] == a[y-1][x] + GAP:
        (allY, allX) = alignUpMost(x, y-1, a, allY, allX)
        allX.append('-')
        allY.append(SEQ2[y-1])
        return (allY, allX)
    elif y > 0 and x > 0 and a[y][x] == a[y-1][x-1] + calcPenalty(x-1, y-1):
        (allY, allX) = alignUpMost(x-1, y-1, a, allY, allX)
        allX.append(SEQ1[x-1])
        allY.append(SEQ2[y-1])
        return (allY, allX)
    else:
        (allY, allX) = alignUpMost(x-1, y, a, allY, allX)
        allX.append(SEQ1[x-1])
        allY.append('-')
        return (allY, allX)

# Backtrack steps for best global downmost alignments
# x, y: current coordinates (recursion)
# a: similarity matrix
# returns tuple of lists (alignments)
def alignDownMost(x, y, a, allY=[], allX=[]):
    if x == 0 and y == 0:
        return (allY, allX)
    elif x > 0 and a[y][x] == a[y][x-1] + GAP:
        (allY, allX) = alignDownMost(x-1, y, a, allY, allX)
        allX.append(SEQ1[x-1])
        allY.append('-')
        return (allY, allX)        
    elif y > 0 and x > 0 and a[y][x] == a[y-1][x-1] + calcPenalty(x-1, y-1):
        (allY, allX) = alignDownMost(x-1, y-1, a, allY, allX)
        allX.append(SEQ1[x-1])
        allY.append(SEQ2[y-1])
        return (allY, allX)
    else:
        (allY, allX) = alignDownMost(x, y-1, a, allY, allX)
        allX.append('-')
        allY.append(SEQ2[y-1])
        return (allY, allX)

if __name__ == "__main__":
    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Needleman-Wunsch Script')
    parser.add_argument('seq_1', type=str, help="String of first sequence")
    parser.add_argument('seq_2', type=str, help="String of second sequence")
    parser.add_argument('--gap', nargs='?', const=-2, default=-2, type=int, help="Cost of gaps, default=-2")
    parser.add_argument('--match', nargs='?', const=1, default=1, type=int, help="Reward of match, default=1")
    parser.add_argument('--mismatch', nargs='?', const=-1, default=-1, type=int, help="Cost of mismatch, default=-1")
    parser.add_argument('--down', action='store_true', default=False, help='Use downmost alignment instead of upmost, default=False')
    main(parser.parse_args())

