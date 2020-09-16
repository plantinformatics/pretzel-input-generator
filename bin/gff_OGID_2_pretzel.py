#!/usr/bin/env python3

from __future__ import print_function
from collections import OrderedDict
import argparse, sys, os, json


def eprint(*args, **kwargs):
  print(*args, file=sys.stderr, **kwargs)

parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
  formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=60),
  description="Convert \"similarity\" features from gff to a Pretzel data set JSON")
parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),
  default=sys.stdin, help='provide input file name or use stdin')
parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'),
  default=sys.stdout, help='provide output file name or use stdout')
parser.add_argument('--make-private', action="store_true",
  default=False, help='"Make output dataset private (is public by default)"')
parser.add_argument('--parent', required=True, help='Parent data set name')
parser.add_argument('--name',  required=True, help='Dataset name')
parser.add_argument('--namespace',  required=True, help='Dataset namespace')
parser.add_argument('--short-name', help='Dataset shortName')
parser.add_argument('--source', help='Source of the dataset')
parser.add_argument('--citation', help='Source of the dataset')


args = parser.parse_args()

annotation = OrderedDict()

annotation['name'] = args.name #"${basename}_${seqType}"
annotation['namespace'] = args.namespace
annotation['parent'] = args.parent
annotation['public'] = not args.make_private
annotation['meta'] = OrderedDict()
if args.short_name is not None:
  annotation['meta']['shortName'] = args.short_name
if args.source is not None:
  annotation['meta']['source'] = args.source
if args.citation is not None:
  annotation['meta']['citation'] = args.citation

scope = OrderedDict()

for line in args.infile:
  toks = line.split('\t')
  chr = toks[0]
  #strip prefix
  if chr.lower().startswith("chr"):
    chr = chr[3:]
  #First time seeing a block?
  if chr not in scope:
    scope[chr] = []
  #Extract feature id
  attributes = OrderedDict(s.split('=') for s in toks[8].split(';')[0:2])
  #Add feature to block
  scope[chr].append(OrderedDict([
    ("name", attributes["OGID"]),
    ("value", [ int(toks[3]), int(toks[4]) ]),
  ]))

annotation['blocks'] = []
for key in scope:
  current = OrderedDict( [("scope", key), ("featureType", "linear"), ("features", OrderedDict())] )
  current['features'] = scope[key]
  annotation['blocks'].append(current)

#sort blocks
annotation['blocks'] = sorted(annotation['blocks'],  key=lambda block : block['scope'])

for block in annotation['blocks']:
  #sort by start position in block
  block['features'] = sorted(block['features'],  key=lambda loc : int(loc['value'][0]))


json.dump(annotation, args.outfile, indent=2)
args.outfile.close()