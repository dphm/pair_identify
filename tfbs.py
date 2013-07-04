class TFBS(object):
    def __init__(self, chromosome, code, kind):
        self.code = code
        self.kind = kind
        self.filepath = "/scratch/dpham4/data/%s/%s.txt" % (chromosome, code)
        self.sites = []
        self.num_sites = 0
    
    def fill(self):
        sites = []
        num_sites = 0
        
        with open(self.filepath) as f_in:
            for l in f_in:
                line = l.replace("\x00", "").split()
                
                try:
                    i = int(line[1])
                    if 0 <= i < size:
                        sites.append(i)
                        num_sites += 1
                except ValueError:
                    continue
        
        if self.kind == "MINOR":
            # remove TFBS input file
            call("rm -f %s" % self.filepath, shell=True)
        
        self.sites = sites
        self.num_sites = num_sites