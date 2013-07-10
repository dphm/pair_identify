from collections import defaultdict

def in_chip(data, site1):
    if data[site1] % 10 == 1:
        return True
    else:
        return False

def study(data, tf1, tf2, max_dist):
    d_TTT = []
    d_FTT = []
    freq = defaultdict(int)
    count = defaultdict(int)
    
    tf1_sites = tf1.sites
    tf2_sites = tf2.sites
    
    tf1_len = tf1.num_sites
    tf2_len = tf2.num_sites
    
    l = 0
    r = 0
    
    l_seek = 0
    r_seek = 0
    
    while l < tf1_len and r < tf2_len:
        finished = l_seek >= tf1_len or r_seek >= tf2_len
        
        if not finished:
            site1 = tf1_sites[l_seek]
            site2 = tf2_sites[r_seek]
            dist = abs(site1 - site2)
        
        if not finished and dist <= max_dist:
            # exclude sites in RepeatMasker regions
            if data[site1] < 10 and data[site2] < 10:
                csv_row = "%i,%i,%i" % (site1, site2, dist)
                freq[dist] += 1
        
                if in_chip(data,site):
                    d_TTT.append(csv_row)
                else:
                    d_FTT.append(csv_row)
            
            if site1 >= site2:
                l_seek += 1
            else:
                r_seek += 1
        else:
            if in_chip(data, site1):
                count["TTF"] += 1
            else:
                count["FTF"] += 1
            
            if site1 >= site2:
                r += 1
                r_seek = r
                l_seek = l
            else:
                l += 1
                l_seek = l
                r_seek = r
    
    count["TTT"] = len(d_TTT)
    count["FTT"] = len(d_FTT)
    
    return (d_TTT, d_FTT, freq, count)