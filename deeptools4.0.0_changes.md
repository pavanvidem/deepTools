# computeMatrix

 - --sortRegions 'no' option no longer exists
 - Sorting ascend / descend no longer has subsorting by position.



# normalization

Exactscaling is no longer an option, it's always performed.

# Todo

 - allow multithreaded bw writing
 - properly divide region work over threads -> region sorting & taking size into account
 - calc for computeMatrix functions -> Struct / Enum
 - filehanlder bed file could all be &str, not clones


# Testing
 
## computeMatrix
### reference-point
 - referencePoint: TSS, center, TES
 - sortRegions: descend, ascend, keep
 - sortUsing: mean, median, max, min, sum, region_length
 - averageTypeBins: mean, median, min, max ,std, sum
