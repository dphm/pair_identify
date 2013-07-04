#!/usr/bin/env python

import sys

MAX_ARRAY_SIZE = 250000000

class ARGS(object):
    def __init__(self, chromosome, TF1_name, TF1_code, TF2_code):
        self.chr = chromosome
        self.TF1_name = TF1_name
        self.TF1_code = TF1_code
        self.TF2_code = TF2_code
    
    def __repr__(self):
        return ("%s %s %s %s" % 
               (self.chr, self.TF1_name, self.TF1_code, self.TF2_code))

def all():
    pass

def run():
    pass

def main(argv=sys.argv[1:]):
    pass

if __name__ == '__main__':
    status = main()
    sys.exit(status)