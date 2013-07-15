def read(filepath):
    data = []
    
    with open(filepath) as f_in:
        for l in f_in:
            line = l.split()
            
            try:
                start = int(line[6])
                finish = int(line[7]) + 1
                
                data.append("%i %i\n" % (start, finish))
            except ValueError:
                continue
    
    return data

def write(filepath, data):
    with open(filepath, "w") as f_out:
        for l in data:
            f_out.write(l)
        
def main(argv=sys.argv[1:]):
    chromosome = argv[0]
    filepath = "/scratch/dpham4/PI/data/%s/rmsk.txt" % chromosome
    
    data = read(filepath)
    write(filepath, data)

if __name__ == '__main__':
    main()