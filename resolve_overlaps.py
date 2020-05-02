#!/usr/bin/env python3.6

import argparse
import itertools
import numpy as np
import sys

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
    parser.add_argument('--summarize', action='store_true',
                        help='don\'t output, only summarize results')
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
        array = np.loadtxt(f, dtype=gff_dtype, delimiter='\t', comments='#')
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


def select_sequences(src):
    by_finish = np.sort(src, order=['end'])
    return 0, []


if __name__ == "__main__":
    args = parse_args()

    print('Loading data... ', end='', flush=True, file=sys.stderr)
    array = load_array(args)
    print('%d rows loaded.' % len(array), file=sys.stderr)

    print('Grouping by seqid, type, strand... ', end='', flush=True, file=sys.stderr)
    views = group_views(array)
    print('%d groups found.' % len(views), file=sys.stderr)

    results = []
    for view in views:
        print('Processing group (%s, %s, %s) with %d items... ' %
              (view[0]['seqid'].decode(), view[0]['type'].decode(),
               view[0]['strand'].decode(), len(view)),
              end='', flush=True, file=sys.stderr)
        score, result = select_sequences(view)
        print('score %d, %d rows.' % (score, len(result)), file=sys.stderr)
        results.append(result)

    if not args.summarize:
        print('##gff-version 3')
        for resultset in results:
            for row in resultset:
                print('\t'.join(str(x) for x in row))
