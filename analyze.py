#!/usr/bin/env python

import sys
from rmsk import RMSK
from chipseq import ChipSeq
from tfbs import TFBS
from os import mkdir, path
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
        
        generate_data(args, chromosome, tf1, tf2)
        
        print "Completed: %s" % args

def setup(args):
    chromosome = args.chr
    tf1_name   = args.tf1_name
    tf1_code   = args.tf1_code
    tf2_code   = args.tf2_code
    
    # input filepaths
    rmsk_in = "%s/data/rmsk.txt.gz" % path
    chip_in = ("/home/mcb/blanchem/wgEncodeRegTfbsClustered/%s_%s.bed" %
              (tf1_code[0:6], tf1_name))
    tf1_in  = "/scratch/blanchem/%s/sites/sites.%s.gz" % (chromosome, tf1_code)
    tf2_in  = "/scratch/blanchem/%s/sites/sites.%s.gz" % (chromosome, tf2_code)
    
    # output filepaths
    rmsk_out = "%s/data/%s/rmsk.txt" % (path, chromosome)
    chip_out = "%s/data/%s/chip_seq_%s.txt" % (path, chromosome, tf1_code[0:6])
    tf1_out  = "%s/data/%s/%s.txt" % (path, chromosome, tf1_code)
    tf2_out  = "%s/data/%s/%s.txt" % (path, chromosome, tf2_code)
    
    # prepare textfiles
    if not path.exists(rmsk_out):
        rmsk_status = call("zcat %s | grep -w %s > %s" %
                          (rmsk_in, chromosome, rmsk_out), shell=True)
        if rmsk_status != 0:
            return "Error: unable to prepare %s" % rmsk_out
    
    if not path.exists(chip_out):
        chip_status = call("grep -w %s %s > %s" %
                          (chromosome, chip_in, chip_out), shell=True)
        if chip_status != 0:
            return "Error: unable to prepare %s" % chip_out
    
    if not path.exists(tf1_out):
        tf1_status = call("zcat %s | grep -w '^0' > %s" %
                         (tf1_in, tf1_out), shell=True)
        if tf1_status != 0:
            return "Error: unable to prepare %s" % tf1_out
    
    if not path.exists(tf2_out):
        tf2_status = call("zcat %s | grep -w '^0' > %s" %
                         (tf2_in, tf2_out), shell=True)
        if tf2_status != 0:
            return "Error: unable to prepare %s" % tf2_out
    
    return "Success"

def generate_data(args, chromosome, tf1, tf2):
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