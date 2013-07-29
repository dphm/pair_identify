pair_identify
=============
Study of transcription factor cooperativity

# Understanding the data

The files used and generated by *analyze.py* can be accessed through the **koksoak.cs.mcgill.ca** server in the **/scratch/dpham4/PI/data** directory.

## Input files

### tf1_list.txt
**Path:** /scratch/dpham4/PI/data/tf1_list.txt  
**Format:** [tf1 code] [tf1 name]\\n  
**Example:** M00008 SP1  

**Description:**  
This file is a list of tf1s (transcription factors with available ChIP-seq data). M\* indicates the forward strand, whereas M\*\_r indicates the reverse strand. This list was compiled from filenames matching M\*_\*.bed in **/home/mcb/blanchem/wgEncodeRegTfbsClustered**. Some tf1s were removed that contained noPhastCons in the filename, or due to empty TFBS files.  

When running *analyze.py*, lines beginning with the \# symbol will be ignored.

### tf2_list.txt

**Path:** /scratch/dpham4/PI/data/[chromosome]/tf2_list.txt  
**Format:** [tf2 code]\\n  
**Example:** M00001  

**Description:**  
This file is a list of tf2s (all transcription factors being tested as candidate partners for tf1s). M\* indicates the forward strand, whereas M\*\_r indicates the reverse strand. This list was compiled from filenames matching sites.M\*.gz in **/scratch/blanchem/[chromsome]/sites**. Some tf2s were removed from the list due to empty TFBS files.

### rmsk.txt : RepeatMasker data

**Path:** /scratch/dpham4/PI/data/[chromosome]/rmsk.txt  
**Format:** [start position],[end position]\\n  
**Example:** 10000 10469  

**Description:**  
This file contains the positions of RepeatMasker regions in a particular chromosome. RepeatMasker data provides regions of known transposable elements in the genome. Extracted from the RepeatMasker data for all chromosomes at **/scratch/dpham4/PI/data/rmsk.txt.gz**.  

### chip_seq_[tf1 code].txt : ChIP-seq data

**Path:** /scratch/dpham4/PI/data/[chromosome]/chip_seq_[tf1 code].txt  
**Format:** [chromosome] [start position] [end position] [tf1 name]\\n  
**Example:** chr1	713863	714256	SP1  

**Description:**  
This file contains the positions of ChIP-seq regions for a particular tf1 in a particular chromosome. These are experimentally verified regions that have been sequenced after being 'fished out' using the given tf1. The tf1 may not be directly bound to the region. Extracted from **/home/mcb/blanchem/wgEncodeRegTfbsClustered/[tf1 code]_[tf1 name].bed**.

### M\*.txt : TFBS data

**Path:** /scratch/dpham4/PI/data/[chromosome]/[tf1 or tf2 code].txt  
**Format:** [genome] [site position] [score] [sequence]\\n  
**Example:** 0	10974	-9.58	GAGGCGTGGC  

**Description:**  
This file contains TFBS positions for either a tf1 or tf2. The 0 in the first column indicates that the sites are found in the human genome (hg19). The site position in the second column is used in *analyze.py*, but the remaining columns are unused. Extracted from **/scratch/blanchem/[chromsome]/sites/sites.[tf1 or tf2 code].gz**

## Output files  

### d_TTT_[tf2 code].csv : TFBS pairs and distances
#### Inside tf1 ChIP-seq regions

**Path:** /scratch/dpham4/PI/data/[chromosome]/[tf1 code]/d_TTT_[tf2 code].csv  
**Format:** [tf1 site position],[tf2 site position],[distance]\\n  
**Example:** 23181478,23181471,7  

**Description:**  
This file reports the positions of pairs with predicted sites that are found in tf1 ChIP-seq regions, and the absolute value between them, if the distance between pairs is less than MAX_TFBS_DIST. The positions may be useful for searching the [UCSC Genome Browser](http://genome.ucsc.edu/cgi-bin/hgGateway) or the [Bejerano GREAT tool](http://bejerano.stanford.edu/great/public/html/).

### d_FTT_[tf2 code].csv : TFBS pairs and distances
#### Outside tf1 ChIP-seq regions

**Path:** /scratch/dpham4/PI/data/[chromosome]/[tf1 code]/d_FTT_[tf2 code].csv  
**Format:** [tf1 site position],[tf2 site position],[distance]\\n  
**Example:** 58313,58217,96  

**Description:**  
This file reports the positions of pairs predicted sites that are not found in tf1 ChIP-seq regions, and the absolute value between them, if the distance between pairs is less than MAX_TFBS_DIST. The positions may be useful for searching the [UCSC Genome Browser](http://genome.ucsc.edu/cgi-bin/hgGateway) or the [Bejerano GREAT tool](http://bejerano.stanford.edu/great/public/html/).

### f_[tf2 code].csv : distance frequencies

**Path:** /scratch/dpham4/PI/data/[chromosome]/[tf1 code]/f_[tf2 code].csv  
**Format:** [distance],[frequency],[z-score]\\n  
**Example:** 7,929,155.999906228  

**Description:**  
This file reports (regardless of TTT or FTT categorization) the distance between all pairs, the number of times that distance was found, and a z-score. The z-score is calculated as (frequency - mean) / (standard deviation), where the mean and standard deviation have been calculated for distances greater than or equal to MIN_MEAN_CUTOFF.  

The purpose of MIN_MEAN_CUTOFF is to determine a 'background level' frequency. The first few distances are ignored as they may have large frequencies from similar PWMs (Position Weight Matrices). See Professor Blanchette for an explanation of the similar PWMs issue.

### s_[tf2 code].txt : number of sites and cases

**Path:** /scratch/dpham4/PI/data/[chromosome]/[tf1 code]/s_[tf2 code].txt  
**Format:**  
\<tf1 code\> \<number of tf1 sites\>\\n  
\<tf2 code\> \<number of tf2 sites\>\\n  
TTT \<number of TTT cases\>\\n  
FTT \<number of FTT cases\>\\n  
**Example:**  
M00795_r 36683  
M00059_r 84099  
TTT 7  
FTT 4349  

**Description:**  
This file reports counts of the number of predicted tf1 and tf2 sites, as well as the number of pairs categorized as TTT or FTT.

### z.csv : pairs with z-scores greater or equal to Z_THRESHOLD

**Path:** /scratch/dpham4/PI/data/z.csv  
**Format:** [args],[highest z-score]\\n  
**Example:** chr1 POU2F2 M00795_r M00059_r,155.999906228  

**Description:**  
A line is appended to this file if the pair's frequency data contains a z-score greater or equal to Z_THRESHOLD. The line contains args (chromosome, tf1 name, tf1 code, tf2 code), and the highest z-score found.

# Usage

*Instructions coming soon!*

---

Feel free to use or modify any part of the existing code to suit the needs of your project.

For questions or clarification, email **danielle.pham@mail.mcgill.ca**.