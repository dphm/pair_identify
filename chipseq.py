class ChipSeq(object):
    def __init__(self, chromosome, code):
        self.code = "chip_seq_" + code[0:6]
        self.filepath = ("/scratch/dpham4/data/%s/%s.txt" %
                        (chromosome, self.code))
    
    def fill(self, data):
        with open(self.filepath) as f_in:
            for l in f_in:
                line = l.split()
                
                try:
                    start = int(line[1])
                    finish = int(line[2]) + 1
                    
                    for i in xrange(start, finish):
                        data[i] += 1
                except ValueError:
                    continue