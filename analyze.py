#!/usr/bin/env python

import os
import sys
import time
from rmsk import RMSK
from chipseq import ChipSeq
from tfbs import TFBS
from casedata import study
from stats import z_scores
from subprocess import call

path = "/scratch/dpham4/PI/data"
log = "log"

MAX_TFBS_DIST = 50
MIN_MEAN_CUTOFF = 20
Z_THRESHOLD = 5

class ARGS(object):
    """Store information about chromosome and pair of tfs being compared"""
    def __init__(self, chromosome, tf1_name, tf1_code, tf2_code):
        self.chr = chromosome
        self.tf1_name = tf1_name
        self.tf1_code = tf1_code
        self.tf2_code = tf2_code
    
    def __repr__(self):
        """Return string representation of args"""
        return ("%s %s %s %s" % 
               (self.chr, self.tf1_name, self.tf1_code, self.tf2_code))
        

def get_time():
    """Return 24h string representation of local time (HH:MM:SS)"""
    return time.strftime("%H:%M:%S", time.localtime())

def tf_lists(chromosome):
    """Return tf lists after reading from files"""
    tf1path = "%s/tf1_list.txt" % path
    tf2path = "%s/%s/tf2_list.txt" % (path, chromosome)
    
    tf1_list = {} # key: code, value: name
    tf2_list = []
    
    with open(tf1path) as list_1:
        # MXXXXX NAME\n
        for entry in list_1:
            if entry[0] == "#":
                continue
            
            code, name = entry.split()
            tf1_list[code] = name.rstrip()
    
    with open(tf2path) as list_2:
        # MXXXXX\n
        for entry in list_2:
            name = entry.rstrip()
            tf2_list.append(name)
    
    return tf1_list, tf2_list

def run(args, rmsk, chip, tf1):
    """Complete data loading and begin data generation"""
    with open(log, "a") as f_log:
        msg = "Processing (%s): %s\n" % (get_time(), args)
        f_log.write(msg)
    
        tf2 = TFBS(args.chr, args.tf2_code)
        generate_data(args, rmsk, chip, tf1, tf2)
    
        msg = "Completed (%s): %s\n" % (get_time(), args)
        f_log.write(msg)

def generate_data(args, rmsk, chip, tf1, tf2):
    """Generate distance, frequency, Z-score data and write to files"""
    chromosome = args.chr
    tf1_code = tf1.code
    tf2_code = tf2.code
    
    # generated data
    d_TTT, d_FTT, freq = study(rmsk, chip, tf1, tf2, MAX_TFBS_DIST)
    z = z_scores(freq, MIN_MEAN_CUTOFF)
    
    # output files
    filepath = "%s/%s/%s" % (path, chromosome, tf1_code)
    sitepath = "%s/s_%s.txt" % (filepath, tf2_code)
    dTTTpath = "%s/d_TTT_%s.csv" % (filepath, tf2_code)
    dFTTpath = "%s/d_FTT_%s.csv" % (filepath, tf2_code)
    freqpath = "%s/f_%s.csv" % (filepath, tf2_code)
    zpath = "%s/z.csv" % path
    
    # create directory for tf1 output files
    if not os.path.exists(filepath):
        os.mkdir(filepath)
    
    # write number of sites and cases
    with open(sitepath, "w") as f_out:
        f_out.write("%s %s\n" % (tf1_code, tf1.num_sites))
        f_out.write("%s %s\n" % (tf2_code, tf2.num_sites))
        f_out.write("%s %i\n" % ("TTT", len(d_TTT)))
        f_out.write("%s %i\n" % ("FTT", len(d_FTT)))
    
    # write positions and distances
    with open(dTTTpath, "w") as f_TTT:
        for d in d_TTT:
            f_TTT.write(d)
    
    with open(dFTTpath, "w") as f_FTT:
        for d in d_FTT:
            f_FTT.write(d)
    
    # write frequencies and Z-scores
    with open(freqpath, "w") as f_freq:
        for d in freq:
            score = None
            if d in z:
                score = z[d]
            
            f_freq.write("%i,%i,%s\n" % (d, freq[d], score))
    
    highscore = 0
    
    if z:
        highscore = max(z.values())
    
    if highscore >= Z_THRESHOLD:
        with open(log, "a") as f_log:
            msg = "%s,%s\n" % (args, highscore)
            f_log.write(msg)
    
def main(argv=sys.argv[1:]):
    # command line processing
    arg_len = len(argv)
    
    # all pairs of tfs
    if arg_len == 2:
        if argv[0] == "--all":
            call("echo -n > log", shell=True) # clear log
            
            chromosome = argv[1]
            
            rmsk = RMSK(chromosome)
            tf1_list, tf2_list = tf_lists(chromosome)
    
            for tf1_code in tf1_list:
                tf1_name = tf1_list[tf1_code]
            
                chip = ChipSeq(chromosome, tf1_code)
                tf1 = TFBS(chromosome, tf1_code)
        
                for tf2_code in tf2_list:
                    if tf1_code != tf2_code:
                        args = ARGS(chromosome, tf1_name, tf1_code, tf2_code)
                        run(args, rmsk, chip, tf1)
            
            return 0
        else:
            print "usage:   python analyze.py --all CHR"
            print "example: python analyze.py --all chr1"
            return 1
    # one pair of tfs
    elif arg_len == 4:
        chromosome = argv[0]
        tf1_name = argv[1]
        tf1_code = argv[2]
        tf2_code = argv[3]
        
        args = ARGS(chromosome, tf1_name, tf1_code, tf2_code)
        rmsk = RMSK(chromosome)
        chip = ChipSeq(chromosome, tf1_code)
        tf1 = TFBS(chromosome, tf1_code)
        
        run(args, rmsk, chip, tf1)
        
        return 0
    else:
        print "usage:   python analyze.py CHR TF1_NAME TF1_CODE TF2_CODE"
        print "example: python analyze.py chr1 NFKB M00774 M00497"
        return 1

if __name__ == '__main__':
    status = main()
    sys.exit(status)