#!/usr/bin/env python3.6

# Runtime on my machine: ~6mins parsing, ~1.5mins processing.
# Without ulimit however, it just seems to run out of memory.

# An implementation quirk due to using linked lists to ameliorate the O(n^2)
# blowup in intermediate sequences - and to maximize data sharing - is that the
# output is reversed.

import argparse
from collections import defaultdict
import heapq
import itertools
import sys

import numpy as np

# Pandas would work better, numpy requires fixed-width string columns
# But I didn't find pandas on merlin...
gff_dtype = np.dtype([('seqid', '|S10'), ('source', '|S15'), ('type', '|S15'),
                      ('start', int), ('end', int), ('score', float),
                      ('strand', 'a1'), ('phase', 'a1')])


def parse_args():
    parser = argparse.ArgumentParser(description='Sequence overlap solver')
    parser.add_argument('input', metavar='INPUT', help='input GFF file')
    parser.add_argument('-s', '--sorted', action='store_true',
                        help='indicates the input is already sorted by start time')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='don\'t output GFF, only summarize results')
    parser.add_argument('--cache', action='store_true',
                        help='use a pickled .npy as a cache')
    return parser.parse_args()


def load_array(args):
    '''Load GFF into a numpy array and optionally cache the parsed data'''
    if args.cache:
        try:
            return np.load(args.input + '.npy')
        except IOError:
            pass
    with open(args.input, 'r') as f:
        array = np.loadtxt(f, dtype=gff_dtype)
        # Sort by start time if not already pre-sorted
        if not args.sorted:
            # In-place sort
            array.sort(order=['start'])
        if args.cache:
            np.save(args.input, array)
        return array


def group_views(src):
    '''Group by (seqid, type, strand) into views'''
    views = []
    xs, ys, zs = np.unique(src['seqid']), np.unique(src['type']), np.unique(src['strand'])
    for a, b, c in itertools.product(xs, ys, zs):
        views.append(src[np.logical_and.reduce([
            src['seqid'] == a, src['type'] == b, src['strand'] == c])])
    return views


def prepare_p_set(by_end):
    '''Prepare the p(i) set/function of closest previous non-overlapping sequence'''
    by_start = np.argsort(by_end, order=['end'])
    # p(i), I guess the closest to a "linear pass" I can get?
    print('Preprocessing p(i)... ', end='', flush=True, file=sys.stderr)
    p = defaultdict(lambda: -1)
    p_inprogress = []
    for seq in reversed(by_start):
        while p_inprogress and by_end[seq]['end'] < -p_inprogress[0][0]:
            _, seqnum = heapq.heappop(p_inprogress)
            p[seqnum] = seq
        heapq.heappush(p_inprogress, (-by_end[seq]['start'], seq))
    print('done.', file=sys.stderr)
    return p

def select_sequences(by_start):
    '''Main body of the algorithm'''
    by_end = np.sort(by_start, order=['end'])
    p = prepare_p_set(by_end)

    # Avoiding the O(n^2) seems to need a recursive (as opposed to an iterative) solution
    print('Calculating optimal sets... ', end='', flush=True, file=sys.stderr)
    a = { -1: [], 0: (0, None) }
    s = { -1: 0, 0: by_end[0]['score'] }
    best_a, best_s = -1, -1
    for i in range(0, len(by_end)):
        w = by_end[i]['score']
        if s[i - 1] > w + s[p[i]]:
            a[i] = a[i - 1]
            s[i] = s[i - 1]
        else:
            a[i] = (i, a[p[i]])  # Not as inefficient due to sharing, I guess
            s[i] = s[p[i]] + w
            if s[i] > best_s:
                best_s, best_a = s[i], a[i]
    print('done.', file=sys.stderr)

    return best_s, best_a, by_end


def main():
    args = parse_args()

    print('Loading data... ', end='', flush=True, file=sys.stderr)
    array = load_array(args)
    print('%d rows loaded.' % len(array), file=sys.stderr)

    print('Grouping by seqid, type, strand... ', end='', flush=True, file=sys.stderr)
    views = group_views(array)
    print('%d groups found.' % len(views), file=sys.stderr)

    if not args.quiet:
        print('##gff-version 3')

    for view in views:
        print('Processing group (%s, %s, %s) with %d items... ' %
              (view[0]['seqid'].decode(), view[0]['type'].decode(),
               view[0]['strand'].decode(), len(view)),
              file=sys.stderr)
        score, result, src = select_sequences(view)
        print('Best score %d, %d rows.' % (score, len(result)), file=sys.stderr)

        if not args.quiet:
            while result:
                row = src[result[0]]
                result = result[1]
                print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
                    row['seqid'].decode(), row['source'].decode(),
                    row['type'].decode(), row['start'], row['end'],
                    row['score'], row['strand'].decode(),
                    row['phase'].decode()))


if __name__ == "__main__":
    main()
