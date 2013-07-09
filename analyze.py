#!/usr/bin/env python

import os
import sys
import time
import multiprocessing as mp
import stats
import casedata
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
    call("echo -n > log", shell=True) # clear log
    all_args = get_args(argv)
    
    mgr = mp.Manager()
    q = mgr.Queue()
    pool = mp.Pool(processes=6, maxtasksperchild=1)
    
    appender = pool.apply_async(file_append, (q,))
    jobs = []
    
    for args in all_args:
        print args
        # job = pool.apply_async(run, (q, args))
        # jobs.append(job)
    
    for job in jobs:
        job.get()
    
    q.put(None)
    pool.close()

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
            f_out.write("%s\n" % message)
        
        time.sleep(0.001)

def run(q, argv):
    args = ARGS(argv[0], argv[1], argv[2], argv[3])
    log = "log"
    
    q.put(log, "Processing: %s" % args)
    
    # file preparation
    status = setup(args)
    
    if status != "Success":
        q.put(log, status)
        q.put(log, "Failed: %s" % args)
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
    
    # fill TFBS lists
    tf1.fill()
    tf2.fill()
    
    generate_data(q, args, args.chr, rmsk, chip, tf1, tf2)
    
    q.put(log, "Completed: %s" % args)

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

def generate_data(q, args, chromosome, rmsk, chip, tf1, tf2):
    d_TTT, d_FTT, freq, count = study(data, tf1, tf1, max_dist)
    
    for key in freq:
        print "freq[%i] = %i" % key, freq[key]
    
    for key in count:
        print "count[%s] = %i" % key, count[key]
    
    # z-scores
    # write to files

def main(argv=sys.argv[1:]):
    # command line processing
    arg_len = len(argv)
    
    if arg_len == 4:
        mgr = mp.Manager()
        q = mgr.Queue()
        pool = mp.Pool(processes=2, maxtasksperchild=1)
        
        appender = pool.apply_async(file_append, (q,))
        
        job = pool.apply_async(run, (q, argv))
        job.get()
        
        q.put(None)
        pool.close()
        
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