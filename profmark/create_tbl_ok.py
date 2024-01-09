#! /usr/bin/env python3

import sys
import argparse
import random

usage = 'create_tbl_ok.py <tblfile> <tblfile_ok>'
#if len(sys.argv) != 3: sys.exit('Incorrect number of cmdline args.\nUsage: {}'.format(usage))
#(tblfile, tblfile_ok) = sys.argv[1:]

parser = argparse.ArgumentParser(prog='create_tbl_ok', description='profmark hmmsearch')
parser.add_argument('tblfile', type=str,
                    help='tblfile from profmark')
parser.add_argument('tblfile_ok', type=str,
                    help='tblfile with ok msas')

# options
parser.add_argument('--getfirst',  type=int, default=0, help='take n random msa from tbl')
parser.add_argument('--random',    type=int, default=0, help='take n random msa from tbl')
parser.add_argument('--name',      type=str, default=None, help='take one particular msa from tbl')

args = parser.parse_args()
if (len(sys.argv) < 3):
    parser.print_help()
    sys.exit('Incorrect number of cmdline args.\nUsage: {}'.format(usage))
tblfile     = args.tblfile
tblfile_ok  = args.tblfile_ok

n_msa_ok = 0
with open(tblfile) as tblfp:
    for line in tblfp:
        if line[0] == '#': continue
        
        fields = line.split()
        if fields[7] == 'ok':
            n_msa_ok += 1

msa_rand = None
if (args.random > 0):
    msa_rand = random.sample(range(1, n_msa_ok), args.random)
    
# Read the master table, for MSAs with successful splits
#
msaname = []
alen    = {}
msainfo = {}
n_ok    = 0
n       = 0
with open(tblfile) as tblfp:
    for line in tblfp:
        if line[0] == '#': continue   
        fields = line.split()
        if fields[7] == 'ok':
            n_ok += 1
        
            if args.name:
                if fields[0] == args.name:
                    msaname.append(fields[0])
                    alen[fields[0]]    = int(fields[2])
                    msainfo[fields[0]] = line.strip('\n\r')
                    break
           
            elif args.getfirst > 0:
                if n < args.getfirst:
                    msaname.append(fields[0])
                    alen[fields[0]]    = int(fields[2])
                    msainfo[fields[0]] = line.strip('\n\r')
                    n += 1
                
            elif (args.random > 0):
                if n_ok in msa_rand:
                    msaname.append(fields[0])
                    alen[fields[0]]    = int(fields[2])
                    msainfo[fields[0]] = line.strip('\n\r')
                    
            else:
                msaname.append(fields[0])
                alen[fields[0]]    = int(fields[2])
                msainfo[fields[0]] = line.strip('\n\r')
            

# Sort the list of MSA names by length of alignment (in columns),
# to help with load balancing; spread largest MSAs across the subtbls
#
msaname.sort(key=lambda s:alen[s], reverse=True)
print("msa", len(msaname))
print("tblfile", tblfile_ok)
with open(tblfile_ok, 'w') as f:
    for i in range(len(msaname)):
        print(msainfo[msaname[i]], file=f)
