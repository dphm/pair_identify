from subprocess import call

class TFBS(object):
    def __init__(self, chromosome, code):
        self.code = code
        self.filepath = "/scratch/dpham4/PI/data/%s/%s.txt" % (chromosome, code)
        self.sites = None
        self.num_sites = 0
        
        self.fill()
    
    def fill(self):
        sites = []
        
        with open(self.filepath) as f_in:
            for l in f_in:
                line = l.replace("\x00", "").split()
                
                try:
                    i = int(line[1])
                    if i >= 0:
                        sites.append(i)
                except (ValueError, IndexError) as e:
                    continue
        
        self.sites = sites
        self.num_sites = len(sites)