#!/usr/bin/env python
from __future__ import print_function
import argparse, sys, os, json, gzip

def eprint(*args, **kwargs):
  print(*args, file=sys.stderr, **kwargs)

parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
  formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=60),
  description="Extract a value for a key from top level of a JSON")
parser.add_argument('-k', '--key', default='name',  help='key to retrieve value for, defaults to "name"')
parser.add_argument('i', nargs='+', type=argparse.FileType('r'),
  default='sys.stdin', metavar='INFILE', help='provide input file name or use stdin')

args = parser.parse_args()

for f in args.i:
  if f.name.endswith('.gz'):
    with gzip.GzipFile(f.name, 'r') as fin:
      jsondict = json.loads(fin.read().decode('utf-8'))
  else:
    jsondict = json.load(f)
  if(args.key in jsondict):
    print(jsondict[args.key])
  else:
    eprint("key: "+args.key+" not in top level of "+args.i)
