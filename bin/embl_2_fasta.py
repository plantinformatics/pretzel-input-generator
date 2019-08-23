#!/usr/bin/env python
from __future__ import print_function
import argparse, sys, os
from Bio import SeqIO

def eprint(*args, **kwargs):
  print(*args, file=sys.stderr, **kwargs)

parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
  formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=60),
# parser = argparse.ArgumentParser(os.path.basename(__file__),
  description="Extract fasta from ENA embl e.g. from the output of wget -O LS480641.1.embl \"http://www.ebi.ac.uk/ena/data/view/LS480641.1&display=text\"")
parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),
  default=sys.stdin, help='provide input file name or use stdin')
parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'),
  default=sys.stdout, help='provide output file name or use stdout')


args = parser.parse_args()

count = SeqIO.convert(args.infile, "embl", args.outfile, "fasta")

eprint(count)