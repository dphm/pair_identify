class RMSK(object):
    def __init__(self, chromosome):
        self.filepath = "/scratch/dpham4/PI/data/%s/rmsk.txt" % chromosome
        self.data = None
        
        self.fill()
    
    def fill(self):
        data = set()
        
        with open(self.filepath) as f_in:
            for l in f_in:
                line = l.split()
                
                try:
                    start = int(line[0])
                    finish = int(line[1])
                
                    for i in xrange(start, finish):
                        data.add(i)
                except ValueError:
                    continue
        
        self.data = data