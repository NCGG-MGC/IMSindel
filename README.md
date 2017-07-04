# IMSIndel

## depend tools

* samtools
* mafft
* glsearch
* BioRuby

## usage

```
  ruby bin/imsindel --bam foo.bam --chr 1 --outd out --indelsize 10000 --reffa ref.fa --glsearch-mat mydna.mat -temp /dev/shm/
```

## options

* --bam /path/to/foo.bam
* --chr chromosome
* --outd /path/to/outoput-dir
* --indelsize indel-size
* --reffa /path/to/ref.fa
* --baseq avg-base-quality
* --mapq mapping-quality
* --within read-grouping-bases
* --pair-within grouping-pair-bases
* --alt-read-depth alt-read-depth
* --support-reads grouping-minimum-reads
* --clip-length clipping-flagment-minimum-bases
* --glsearch /path/to/glsearch36
* --glsearch-mat /path/to/foo.mat
* --mafft /path/to/mafft
* --samtools /path/to/samtools
* --temp temp-dir-for-mafft-and-glsearch
* --thread number-of-mafft-threads
