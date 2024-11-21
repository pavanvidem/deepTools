# Counting

Counting reads has slightly changed.  
for -bs 1, the counting mechanism remains the same.  
for -bs > 1, a read is split in multiple contiguous blocks  
  - multiple blocks in 1 bin only count as 1  
  - multiple blocks in multiple bins count as 1 per bin  
  - one block spanning multiple bins counts as 1 in each bin  

# normalization

Exactscaling is no longer an option, it's always performed.

# Todo

 - allow multithreaded bw writing
 - properly divide region work over threads -> region sorting & taking size into account
 