class RMSK(object):
    def __init__(self, chromosome):
        self.filepath = "/scratch/dpham4/PI/data/%s/rmsk.txt" % chromosome
    
    def fill(self, data):
        with open(self.filepath) as f_in:
            for l in f_in:
                line = l.split()
                
                try:
                    start = int(line[6])
                    finish = int(line[7]) + 1
                
                    for i in xrange(start, finish):
                        data[i] += 10
                except ValueError:
                    continue