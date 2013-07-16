#!/usr/bin/env python

import os
import sys
import time
import multiprocessing as mp
from rmsk import RMSK
from chipseq import ChipSeq
from tfbs import TFBS
from casedata import study
from stats import z_scores
from subprocess import call

path = "/scratch/dpham4/PI/data"

MAX_TFBS_DIST = 50
MIN_MEAN_CUTOFF = 20
Z_THRESHOLD = 5

class ARGS(object):
    def __init__(self, chromosome, tf1_name, tf1_code, tf2_code):
        self.chr = chromosome
        self.tf1_name = tf1_name
        self.tf1_code = tf1_code
        self.tf2_code = tf2_code
    
    def __repr__(self):
        return ("%s %s %s %s" % 
               (self.chr, self.tf1_name, self.tf1_code, self.tf2_code))


class DATA(object):
    def __init__(self, chromosome):
        self.rmsk = RMSK(chromosome)
        self.chip = None
        self.tf1 = None
        self.tf2 = None
        
    def load_chip(self, chromosome, tf1_code):
        self.chip = ChipSeq(chromosome, tf1_code)
    
    def load_tf1(self, chromosome, tf1_code):
        self.tf1 = TFBS(chromosome, tf1_code, 1)
    
    def load_tf2(self, chromosome, tf2_code):
        self.tf2 = TFBS(chromosome, tf2_code, 2)


class MP(object):
    def __init__(self, ps):
        self.mgr = mp.Manager()
        self.q = self.mgr.Queue()
        self.pool = mp.Pool(processes=ps, maxtasksperchild=1)
    def activate(self, all_args, data):
        appender = self.pool.apply_async(file_append, (self.q,))
        jobs = []
        
        for args in all_args:
            if data.tf1 == None or data.tf1.code != args.tf1_code:
                data.load_chip(args.chr, args.tf1_code)
                data.load_tf1(args.chr, args.tf1_code)
            
            job = self.pool.apply_async(run, (self.q, args, data))
            jobs.append(job)
        
        for job in jobs:
            job.get()
            
        self.q.put(None)
        self.pool.close()


def all(argv, data):
    call("echo -n > log", shell=True) # clear log
    all_args = get_args(argv)
    
    pool = MP(6)
    pool.activate(all_args, data)

def get_args(argv):
    all_args = []
    
    chromosome = argv[1]
    chippath = "%s/chip_seq_list.txt" % path
    tfbspath = "%s/%s/tfbs_list.txt" % (path, chromosome)
    
    with open(chippath) as chip_list:
        # MXXXXX_NAME.bed
        for chip in chip_list:
            tf1_code, tf1_name = chip.split()

            with open(tfbspath) as tfbs_list:
                # sites.MXXXXX.gz
                for tfbs in tfbs_list:
                    tf2_code = tfbs
                    
                    if tf1_code != tf2_code:
                        all_args.append(ARGS(chromosome, tf1_name,
                                             tf1_code, tf2_code))
    
    return all_args

def file_append(q):
    while 1:
        data = q.get()
        
        if data == None:
            break
        
        file, message = data
        
        with open(file, "a") as f_out:
            f_out.write(message)
        
        time.sleep(0.001)

def run(q, args, data):
    log = "log"
    
    q.put((log, "Processing: %s\n" % args))
    
    # create tf2 object
    data.load_tf2(args.chr, args.tf2_code)
    
    generate_data(q, args, data)
    
    q.put((log, "Completed: %s\n" % args))

def generate_data(q, args, data):
    chromosome = args.chr
    tf1_code = args.tf1_code
    tf2_code = args.tf2_code
    
    # generated data
    d_TTT, d_FTT, freq = study(data, MAX_TFBS_DIST)
    z = z_scores(freq, MIN_MEAN_CUTOFF)
    
    # output files
    filepath = "%s/%s/%s" % (path, chromosome, tf1_code)
    sitepath = "%s/s_%s.txt" % (filepath, tf2_code)
    dTTTpath = "%s/d_TTT_%s.csv" % (filepath, tf2_code)
    dFTTpath = "%s/d_FTT_%s.csv" % (filepath, tf2_code)
    freqpath = "%s/f_%s.csv" % (filepath, tf2_code)
    zpath = "%s/z.txt" % path
    
    if not os.path.exists(filepath):
        os.mkdir(filepath)
    
    # write number of sites and cases
    with open(sitepath, "w") as f_out:
        f_out.write("%s %s\n" % (tf1_code, data.tf1.num_sites))
        f_out.write("%s %s\n" % (tf2_code, data.tf2.num_sites))
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
    
    highscore = max(z.values())
    
    if highscore >= Z_THRESHOLD:
        q.put((zpath, "%s, max Z-score: %s\n" % (args, highscore)))
    
def main(argv=sys.argv[1:]):
    # command line processing
    arg_len = len(argv)
    
    if arg_len == 4:
        args = ARGS(argv[0], argv[1], argv[2], argv[3])
        data = DATA(argv[0])
        
        pool = MP(2)
        pool.activate((args,), data)
        
        return 0
    elif arg_len == 2:
        if argv[0] == "--all":
            data = DATA(argv[1])
            
            all(argv, data)
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