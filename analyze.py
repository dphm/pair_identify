#!/usr/bin/env python

import os
import sys
from rmsk import RMSK
from chipseq import ChipSeq
from tfbs import TFBS
from subprocess import call

path = "/scratch/dpham4/PI"

MAX_ARRAY_SIZE = 250000000

class ARGS(object):
    def __init__(self, chromosome, tf1_name, tf1_code, tf2_code):
        self.chr = chromosome
        self.tf1_name = tf1_name
        self.tf1_code = tf1_code
        self.tf2_code = tf2_code
    
    def __repr__(self):
        return ("%s %s %s %s" % 
               (self.chr, self.tf1_name, self.tf1_code, self.tf2_code))

def all(argv):
    print "all()"

def run(argv):
    args = ARGS(argv[0], argv[1], argv[2], argv[3])
    
    with open("log", "a") as sys.stdout:
        print "Processing: %s" % args
        
        # file preparation
        status = setup(args)
        
        if status != "Success":
            print status
            print "Failed: %s" % args
            return
        
        # create data array
        data = [0] * MAX_ARRAY_SIZE
        
        # create objects
        rmsk = RMSK(args.chr)
        chip = ChipSeq(args.chr, args.tf1_code)
        tf1  = TFBS(args.chr, args.tf1_code, 1)
        tf2  = TFBS(args.chr, args.tf2_code, 2)
        
        # fill data array
        rmsk.fill(data)    # RepeatMasker + 10
        chip.fill(data) # ChIP-Seq     +  1
        
        # fill TFBS lists
        tf1.fill()
        tf2.fill()
        
        generate_data(args, args.chr, rmsk, chip, tf1, tf2)
        
        print "Completed: %s" % args

def setup(args):
    exit_codes = {0: success, 1: rmsk_err, 2: chip_err, 3: tf1_err, 4: tf2_err}
    
    def success():
        return "Success"
    
    def rmsk_err():
        return "Error: unable to prepare %s" % rmsk_out
    
    def chip_err():
        return "Error: unable to prepare %s" % chip_out
    
    def tf1_err():
        return "Error: unable to prepare %s" % tf1_out
    
    def tf2_err():
        return "Error: unable to prepare %s" % tf2_out
        
    status = call("./setup.sh %s" % args, shell=True)
    exit_codes[status]()

def generate_data(args, chromosome, rmsk, chip, tf1, tf2):
    print "generate_data()"

def main(argv=sys.argv[1:]):
    # command line processing
    arg_len = len(argv)
    
    if arg_len == 4:
        run(argv)
        return 0
    elif arg_len == 2:
        if argv[0] == "--all":
            all(argv)
            return 0
        else:
            print "usage:   python analyze.py --all CHR"
            print "example: python analyze.py --all chr1"
            return 1
    else:
        print "usage:   python analyze.py CHR TF1_NAME TF1_CODE TF2_CODE"
        print "example: python analyze.py chr1 NFKB M00774 M00497"
        return 1

if __name__ == '__main__':
    status = main()
    sys.exit(status)