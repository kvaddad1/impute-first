# varcount

## dependencies: 

* htslib
* c++11 compliant compiler.

## Compilation

`make`

## Description and Usage: 

```
./varcount [options] <vcf> <sam>

<vcf>=STR               [bv]cf file name (required)
<sam>=STR               [bs]sam file name (required)
-s/--sample-name=STR    sample name in VCF output
                        (default: sample)
-g/--genotype           [likelihood, alt_sensitive, threshold]. 'predict' a genotype in the GT field
                        likelihood: use crude gt likelihood from counts
                        alt_sensitive: automatically call alleles with any alt evidence as alt/alt
                        threshold: use a manual threshold of the difference between alt and ref count to determine gt. Use -c parameter 
to specifiy threshoold
-c                      int>=0 (default: 0). for use with '-g threshold'. estabilish threshold of ref-alt for determing genotype.  
-a/--min-alt-count      int>=0 (default: 0). filter loci by minimum depth of reads covering alt allele.
-r/--min-ref-count      int>=0 (default: 0). filter loci by minimum depth of reads covering ref allele.
-m/--min-total-count    int>=0 (default: 0). filter loci by minimum depth of total reads.
-k/--keep               ignore filters, print all records regardless of coverage/genotype
-v/--verbose            prints detailed logging information to stderr
-h/--help               print this help message
```
