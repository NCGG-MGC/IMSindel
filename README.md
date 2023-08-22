# IMSIndel

Insertions and deletions (indels) have been implicated in dozens of human diseases through the radical alteration of gene function by short frameshift indels as well as long indels. However, the accurate detection of these indels from next-generation sequencing data is still challenging. This is particularly true for intermediate-size indels (≥50 bp), due to the short DNA sequencing reads. Here, we developed a new method that predicts intermediate-size indels using BWA soft-clipped fragments (unmatched fragments in partially mapped reads) and unmapped reads. 

[Fig1.Shigemizu.pdf](https://github.com/NCGG-MGC/IMSindel/files/12403882/Fig1.Shigemizu.pdf)

## Reference

Daichi Shigemizu, Fuyuki Miya, Shintaro Akiyama, Shujiro Okuda, Keith A Boroevich, Akihiro Fujimoto, Hidewaki Nakagawa, Kouichi Ozaki, Shumpei Niida, Yonehiro Kanemura, Nobuhiko Okamoto, Shinji Saitoh, Mitsuhiro Kato, Mami Yamasaki, Tatsuo Matsunaga, Hideki Mutai, Kenjiro Kosaki & Tatsuhiko Tsunoda,  IMSindel: An accurate intermediate-size indel detection tool incorporating de novo assembly and gapped global-local alignment with split read analysis, Scientific Reports volume 8, Article number: 5608 (2018)
    
## depend tools

* [ruby (2.3 or higher)](https://www.ruby-lang.org/en/)
* [samtools](https://github.com/samtools/samtools)
* [mafft](http://mafft.cbrc.jp/alignment/software/)
* [glsearch](http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml)
* [BioRuby](http://bioruby.org/)

## usage

```console
$ bin/imsindel --bam foo.bam --chr 1 --outd out --indelsize 10000 --reffa ref.fa
```

### run on docker

#### build image

```console
$ git clone https://github.com/NCGG-MGC/IMSindel.git
$ cd IMSindel
$ docker build -t imsindel .
```

#### run imsindel

```console
$ mkdir /path/to/data
$ mv /path/to/your.bam /path/to/data/
$ samtools index /path/to/data/your.bam
$ mv /path/to/ref.fa /path/to/data/
$ docker run --rm -v /path/to/data:/data imsindel --bam /data/your.bam --chr 1 --outd /data --indelsize 10000 --reffa /data/ref.fa
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
* --output-consensus-seq /path/to/output-dir
* --exclude-region /path/to/exclude-list(bed like format 1-started)

## output

#|column|description
-----|------|-----------
1|indel_type|DEL=deletion, INS=insertion
2|call_type|Hete=heterozygous indel(0.15<#indel_depth/#ttl_depth<=0.7), Homo=homozygous indel(#indel_depth/#ttl_depth>0.7)
3|chr|chromosome number
4|sttpos|indel’s start position
5|endpos|indel’s end position
6|indel_length|indel size
7|indel_str|indel sequence
8|#indel_depth|read count including indels
9|#ttl_depth|total read count
10|details(indelcall_indeltype_depth)|composed of four components; <br> 1. Indel_type <br> 2. LI=long insertion, ULI=uncomplete long insertion, LD=long deletion, B: clipped fragments on the right side of read sequences, F: clipped fragments on the left side of read sequences, SI: short indel, <br> 3. #indel_depth, <br> 4. clip_sttpos
11|clip_sttpos|clipped fragments’ start position
12|depth(>=10)|High if #total depth >=10


## update

* 2019/12/19
  * fix fix negative sttpos #10
  * fix fasta url
* 2018/6/13
  * add exclude-region option
  * add mafft error log
