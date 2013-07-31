import sys
import glob

def write(file, tf2_list):
    """Write tf2_list to file"""
    with open(file, 'w') as f_out:
        for tf in tf2_list:
            f_out.write('%s\n' % tf)

def create_tf2_list(chromosome):
    tfbs_path = '/scratch/blanchem/%s/sites' % chromosome
    tf2_path = '/scratch/dpham4/PI/data/%s/tf2_list.txt' % chromosome
    
    tf2_list = []
    
    files = glob.glob('%s/sites.*.gz' % tfbs_path)
    for file in files:
        tf2_code = file.split('/')[5].split('.')[1]
        tf2_list.append(tf2_code)
        
    write(tf2_path, tf2_list)