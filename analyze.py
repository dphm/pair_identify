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
    rmsk = None
    
    def __init__(self):
        self.chip = None
        self.tf1 = None
        self.tf2 = None


class MP(object):
    def __init__(self, ps):
        self.mgr = mp.Manager()
        self.q = self.mgr.Queue()
        self.pool = mp.Pool(processes=ps, maxtasksperchild=1)
    def activate(self, argv):
        log = "log"
        
        appender = self.pool.apply_async(file_append, (self.q,))
        
        if argv[0] == "--all":
            chromosome = argv[1]
            
            msg = "Loading RepeatMasker data (%s)\n" % get_time()
            self.q.put((log, msg))
            
            DATA.rmsk = RMSK(chromosome)
            
            tf1_list, tf2_list = tf_lists(chromosome)
            
            for tf1_code in tf1_list:
                jobs = []
                tf1_name = tf1_list[tf1_code]
                
                msg = "Loading %s data (%s)\n" % (tf1_name, get_time())
                self.q.put((log, msg))
                
                chip = ChipSeq(chromosome, tf1_code)
                tf1 = TFBS(chromosome, tf1_code)
                
                for tf2_code in tf2_list:
                    if tf1_code != tf2_code:
                        args = ARGS(chromosome, tf1_name, tf1_code, tf2_code)
                        
                        data = DATA()
                        data.chip = chip
                        data.tf1 = tf1
                        
                        job = self.pool.apply_async(run, (self.q, args, data))
                        jobs.append(job)
                
                for job in jobs:
                    job.get()
        else:
            chromosome = argv[0]
            tf1_name = argv[1]
            tf1_code = argv[2]
            tf2_code = argv[3]
            
            args = ARGS(chromosome, tf1_name, tf1_code, tf2_code)
            
            DATA.rmsk = RMSK(chromosome)
            
            data = DATA()
            data.chip = ChipSeq(chromosome, tf1_code)
            data.tf1 = TFBS(chromosome, tf1_code)
            
            job = self.pool.apply_async(run, (self.q, args, data))
            job.get()
            
        self.q.put(None)
        self.pool.close()


def get_time():
    return time.strftime("%H:%M:%S", time.localtime())

def file_append(q):
    while 1:
        data = q.get()
    
        if data == None:
            break
    
        file, message = data
    
        with open(file, "a") as f_out:
            f_out.write(message)
    
        time.sleep(0.001)

def tf_lists(chromosome):
    tf1path = "%s/tf1_list.txt" % path
    tf2path = "%s/%s/tf2_list.txt" % (path, chromosome)
    
    tf1_list = {}
    tf2_list = []
    
    with open(tf1path) as list_1:
        # MXXXXX NAME\n
        for entry in list_1:
            code, name = entry.split()
            tf1_list[code] = name.rstrip()
    
    with open(tf2path) as list_2:
        # MXXXXX\n
        for entry in list_2:
            name = entry.rstrip()
            tf2_list.append(name)
    
    return tf1_list, tf2_list

def run(q, args, data):
    log = "log"
    
    msg = "Processing (%s): %s\n" % (get_time(), args)
    q.put((log, msg))
    
    data.tf2 = TFBS(args.chr, tf2_code)
    generate_data(q, args, data)
    
    msg = "Completed (%s): %s\n" % (get_time(), args)
    q.put((log, msg))

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
    
    highscore = 0
    
    if z:
        highscore = max(z.values())
    
    if highscore >= Z_THRESHOLD:
        msg = "%s, max Z-score: %s\n" % (args, highscore)
        q.put((zpath, msg))
    
def main(argv=sys.argv[1:]):
    # command line processing
    arg_len = len(argv)
    
    if arg_len == 2:
        if argv[0] == "--all":
            call("echo -n > log", shell=True) # clear log
            
            pool = MP(9)
            pool.activate(argv)
            
            return 0
        else:
            print "usage:   python analyze.py --all CHR"
            print "example: python analyze.py --all chr1"
            return 1
    elif arg_len == 4:
        pool = MP(2)
        pool.activate(argv)
        
        return 0
    else:
        print "usage:   python analyze.py CHR TF1_NAME TF1_CODE TF2_CODE"
        print "example: python analyze.py chr1 NFKB M00774 M00497"
        return 1

if __name__ == '__main__':
    status = main()
    sys.exit(status)