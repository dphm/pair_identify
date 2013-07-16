class ChipSeq(object):
    def __init__(self, chromosome, code):
        self.code = "chip_seq_" + code[0:6]
        self.filepath = ("/scratch/dpham4/PI/data/%s/%s.txt" %
                        (chromosome, self.code))
        self.sites = None
        
        self.fill()
    
    def fill(self):
        sites = set()
        
        with open(self.filepath) as f_in:
            for l in f_in:
                line = l.split()
                
                try:
                    start = int(line[1])
                    finish = int(line[2]) + 1
                    
                    for i in xrange(start, finish):
                        sites.add(i)
                except ValueError:
                    continue
        
        self.sites = sites