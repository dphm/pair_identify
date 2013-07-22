from collections import defaultdict

def study(rmsk, chip, tf1, tf2, max_dist):
    """Return lists of positions and distances for TTT and FTT cases,
       and dictionary of distance frequencies"""
    d_TTT = []
    d_FTT = []
    freq = defaultdict(int) # key: distance, value: count
    
    rmsk_sites = rmsk.sites
    chip_sites = chip.sites
    tf1_sites = tf1.sites
    tf2_sites = tf2.sites
    
    tf1_len = tf1.num_sites
    tf2_len = tf2.num_sites
    
    l = 0 # bookmark index for tf1_sites
    r = 0 # bookmark index for tf2_sites
    
    l_seek = 0 # comparison index for tf1_sites
    r_seek = 0 # comparison index for tf2_sites
    
    while l < tf1_len and r < tf2_len:
        # no more entries to advance seek index in tf1_sites or tf2_sites
        finished = l_seek >= tf1_len or r_seek >= tf2_len
        
        if not finished:
            site1 = tf1_sites[l_seek]
            site2 = tf2_sites[r_seek]
            dist = abs(site1 - site2)
        
        if not finished and dist <= max_dist:
            # exclude sites in RepeatMasker regions
            if site1 not in rmsk_sites and site2 not in rmsk_sites:
                csv_row = "%i,%i,%i\n" % (site1, site2, dist)
                freq[dist] += 1
        
                if site1 in chip_sites:
                    d_TTT.append(csv_row)
                else:
                    d_FTT.append(csv_row)
            
            if site1 >= site2:
                l_seek += 1
            else:
                r_seek += 1
        else:
            if site1 >= site2:
                r += 1
                r_seek = r
                l_seek = l
            else:
                l += 1
                l_seek = l
                r_seek = r
    
    return (d_TTT, d_FTT, freq)