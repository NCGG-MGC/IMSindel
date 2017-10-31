# IMSIndel

## depend tools

* [samtools](https://github.com/samtools/samtools)
* [mafft](http://mafft.cbrc.jp/alignment/software/)
* [glsearch](http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml)
* [BioRuby](http://bioruby.org/)

## usage

```
  bin/imsindel --bam foo.bam --chr 1 --outd out --indelsize 10000 --reffa ref.fa
```

## options

* --bam /path/to/foo.bam
* --chr chromosome
* --outd /path/to/outoput-dir
* --indelsize maximal-indel-size
* --reffa /path/to/ref.fa
* --baseq avg-base-quality [20]
* --mapq mapping-quality [20]
* --within read-grouping-bases [3]
* --pair-within grouping-pair-bases [5]
* --alt-read-depth alt-read-depth [5]
* --support-reads grouping-minimum-reads [3]
* --clip-length clipped-fragment-minimum-length [5]
* --support-clip-length one-support-length-in-clipped-fragments [5]
* --glsearch /path/to/glsearch36 [glsearch36]
* --glsearch-mat /path/to/foo.mat [data/mydna.mat]
* --mafft /path/to/mafft [mafft]
* --samtools /path/to/samtools [samtools]
* --temp temp-dir-for-mafft-and-glsearch [/tmp]
* --thread number-of-mafft-threads [1]
