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

path = "/scratch/dpham4/PI"

MAX_ARRAY_SIZE = 250000000
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


class MP(object):
    def __init__(self, ps):
        self.mgr = mp.Manager()
        self.q = self.mgr.Queue()
        self.pool = mp.Pool(processes=ps, maxtasksperchild=1)
    def activate(self, all_args):
        appender = self.pool.apply_async(file_append, (self.q,))
        jobs = []
        
        for args in all_args:
            job = self.pool.apply_async(run, (self.q, args))
            jobs.append(job)
        
        for job in jobs:
            job.get()
            
        self.q.put(None)
        self.pool.close()


def all(argv):
    call("echo -n > log", shell=True) # clear log
    all_args = get_args(argv)
    
    pool = MP(6)
    pool.activate(all_args)

def get_args(argv):
    all_args = []
    
    chromosome = argv[1]
    # chippath = "/scratch/dpham4/data/chip_seq_list.txt"
    tfbspath = "/scratch/dpham4/data/%s/tfbs_list.txt" % chromosome
    
    # with open(chippath) as chip_list:
    #    # MXXXXX_NAME.bed
    #    for chip in chip_list:
    #        tf1 = chip.split("_")
    #        tf1_code = tf1[0]
    #        tf1_name = tf1[1].split(".")[0]
    
    tf1_code = "M00975"
    tf1_name = "RFX5"

    with open(tfbspath) as tfbs_list:
        # sites.MXXXXX.gz
        for tfbs in tfbs_list:
            tf2_code = tfbs.split(".")[1]
        
            if tf1_code != tf2_code:
                all_args.append(ARGS(chromosome, tf1_name, tf1_code, tf2_code))
    
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

def run(q, args):
    log = "log"
    
    q.put((log, "Processing: %s\n" % args))
    
    # file preparation
    status = setup(args)
    
    if status != "Success":
        q.put((log, "%s\n" % status))
        q.put((log, "Failed: %s\n" % args))
        return
    
    # create data array
    data = [0] * MAX_ARRAY_SIZE
    
    # create objects
    rmsk = RMSK(args.chr)
    chip = ChipSeq(args.chr, args.tf1_code)
    tf1  = TFBS(args.chr, args.tf1_code, 1)
    tf2  = TFBS(args.chr, args.tf2_code, 2)
    
    # fill data array
    rmsk.fill(data) # RepeatMasker + 10
    chip.fill(data) # ChIP-Seq     +  1
    
    generate_data(q, args, args.chr, data, tf1, tf2)
    
    q.put((log, "Completed: %s\n" % args))

def setup(args):
    def success():
        return "Success"
    
    def rmsk_err():
        return "Error: unable to prepare rmsk file"
    
    def chip_err():
        return "Error: unable to prepare chip file"
    
    def tf1_err():
        return "Error: unable to prepare tf1 file"
    
    def tf2_err():
        return "Error: unable to prepare tf2 file"
        
    exit_codes = {0: success, 1: rmsk_err, 2: chip_err, 3: tf1_err, 4: tf2_err}
        
    status = call("./setup.sh %s" % args, shell=True)
    return exit_codes[status]()

def generate_data(q, args, chromosome, data, tf1, tf2):
    # generated data
    d_TTT, d_FTT, freq = study(data, tf1, tf2, MAX_TFBS_DIST)
    z = z_scores(freq, MIN_MEAN_CUTOFF)
    
    # output files
    filepath = "%s/data/%s/%s" % (path, chromosome, tf1.code)
    sitepath = "%s/s_%s.txt" % (filepath, tf2.code)
    dTTTpath = "%s/d_TTT_%s.csv" % (filepath, tf2.code)
    dFTTpath = "%s/d_FTT_%s.csv" % (filepath, tf2.code)
    freqpath = "%s/f_%s.csv" % (filepath, tf2.code)
    zpath = "%s/data/z.txt" % path
    
    if not os.path.exists(filepath):
        os.mkdir(filepath)
    
    # write number of sites and cases
    with open(sitepath, "w") as f_out:
        f_out.write("%s %s\n" % (tf1.code, tf1.num_sites))
        f_out.write("%s %s\n" % (tf2.code, tf2.num_sites))
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
        
        pool = MP(2)
        pool.activate((args,))
        
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