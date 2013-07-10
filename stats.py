def z_scores(freq, cutoff):
    """Return z_scores dictionary"""
    # number of frequencies considered
    n = sum(1 for d in freq if d >= cutoff)
    if n == 0:
        return {}
    
    # compute mean
    m = mean(freq, n, cutoff)
    
    # compute standard deviation
    sd = std(freq, m, n, cutoff)
    
    # compute and return Z-scores
    z = {} # key = distance between MAJOR and MINOR
           # value = Z-score
    
    for d in freq:
        if sd != 0:
            z[d] = (freq[d] - m) / sd
        else:
            z[d] = None
    
    return z

def mean(freq, n, cutoff):
    """Compute mean"""
    m = sum(float(freq[d]) for d in freq if d >= cutoff)
    return m / n

def std(freq, m, n, cutoff):
    """Compute standard deviation"""
    var = sum((freq[d] - m) ** 2 for d in freq if d >= cutoff)
    return (var / n) ** 0.5