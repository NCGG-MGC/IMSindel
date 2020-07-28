require_relative 'ims_indel/read_collector'
require_relative 'ims_indel/read_group'
require_relative 'ims_indel/consensus'
require_relative 'ims_indel/indel_detector'
require 'pathname'

module IMSIndel
  class App
    def initialize(config)
      @config = default_config.merge(config)
    end

    def check_input
      unless File.exist?(@config[:reffa])
        puts("#{@config[:reffa]}: No such file")
        return false
      end
      unless File.exist?(@config[:bam])
        puts("#{@config[:bam]}: No such file")
        return false
      end
      puts "mafft version:"
      unless system(@config[:mafft], "--version")
        puts("#{@config[:mafft]}: No such file")
        return false
      end
      puts ""
      puts "samtools version:"
      unless system(@config[:samtools], "--version")
        puts("#{@config[:samtools]}: No such file")
        return false
      end
      puts ""
      puts "glsearch version:"
      unless system(@config[:glsearch], "-h")
        puts("#{@config[:glsearch]}: No such file")
        return false
      end
      puts ""
      return true
    end

    def run
      check_input || (return -1)
      reads = ReadCollector.new(@config[:samtools])
      run_collect_indel_reads(reads) || (return 0)
      run_collect_unmapped_reads(reads)
      read_group = run_make_read_groups(reads)
      consensus = run_make_consensus_support_reads(read_group)
      read_group = nil
      run_make_consensus(reads, consensus)
      run_detect_indels(consensus)
      return 0
    end

    def run_collect_indel_reads(reads)
      puts ">1. collecting indel related reads..."

      is_empty = reads.collect_clip_indel_reads(@config[:bam],
                                                @config[:chr],
                                                @config[:mapq],
                                                @config[:baseq],
                                                @config[:clip_length],
                                                @config[:exclude_region])

      puts "#backward_clips:\t#{reads.backward_clips.size}"
      #puts "#pos2max_back_soft_clip:\t#{indel_reads.pos2max_back_soft_clip.size}"
      puts "#forward_clips:\t#{reads.forward_clips.size}"
      #puts "#pos2max_forward_soft_clip:\t#{indel_reads.pos2max_forward_soft_clip.size}"
      puts "#non_clips:\t#{reads.non_clips.size}"
      puts ">1. collecting indel related reads...done\n"
      if is_empty
        puts 'missing indels'
        return false
      end
      return true
    end

    def run_collect_unmapped_reads(reads)
      puts ">2. collecting unmapped reads..."

      reads.collect_unmapped_reads(@config[:bam], @config[:chr], @config[:mapq], @config[:baseq])

      puts "Insert size"
      puts "Avg:\t#{reads.avg_insert_size}"
      puts "SD:\t#{reads.sd}"
      puts "#unmapped reads:\t#{reads.unmapped_reads.size}"
      puts ">2. collecting unmapped reads...done\n"
    end

    def run_make_read_groups(reads)
      puts ">3. considering support reads..."
      read_group = ReadGroup.new(reads, @config[:within], @config[:support_reads], @config[:support_clip_length])

      puts "#backward clip with support reads:\t#{read_group.backward_clips.size}"
      puts "#forward clip with support reads:\t#{read_group.forward_clips.size}"
      puts "#non_clips with suport reads:\t#{read_group.non_clips.size}"
      puts ">3. considering support reads...done\n"
      return read_group
    end

    def run_make_consensus_support_reads(read_group)
      puts ">4. making consensus seqs from support reads..."
      consensus = Consensus.new(@config[:temp], @config[:thread], @config[:mafft], !@config[:output_consensus_seq].nil?)
      consensus.make_support_read_consensus(read_group, @config[:alt_read_depth])
      puts "#backward clip with consensus:\t#{read_group.backward_clips.size} --> #{consensus.backward_clip_consensus.size}"
      puts "#forward clip with consensus:\t#{read_group.forward_clips.size} --> #{consensus.forward_clip_consensus.size}"
      puts "#shot indel with consensus:\t#{read_group.non_clips.size} --> #{consensus.short_indel_consensus.size}"
      puts ">4. making consensus seqs from support reads...done\n"
      return consensus
    end

    def run_make_consensus(reads, consensus)
      puts ">5. making consensus seq from B and F.."
      consensus.make_consensus_seq(reads, @config[:pair_within], @config[:alt_read_depth], @config[:indelsize])
      puts "making consensus seq for long deletion...done\n"
    end

    def run_detect_indels(consensus)
      puts ">6. detection of indels..."
      puts "#paired long indel candidates:\t#{consensus.paired_long_insertions.size}"
      puts "#unpaired long indel candidates:\t#{consensus.unpaired_long_indels.size}"
      puts "#short indel candidates:\t#{consensus.short_indel_consensus.size}"
      outf = Pathname.new(@config[:outd]) + "#{@config[:chr]}.out"
      indel_detector = IndelDetector.new(@config[:glsearch], @config[:glsearch_mat], @config[:samtools], @config[:temp])
      indel_detector.indel_call(consensus.indel_list,
                                @config[:chr],
                                @config[:indelsize],
                                @config[:reffa],
                                @config[:bam],
                                @config[:mapq],
                                @config[:alt_read_depth],
                                outf.to_s)
      puts ">6. detection of indels...done"
      if @config[:output_consensus_seq]
        puts ">7. writing consensus seq..."
        indel_detector.write_indel_reads(@config[:output_consensus_seq], @config[:chr], consensus)
        puts ">7. writing consensus seq...done"
      end
    end

    def print_config
      puts ">Parameters:"
      puts "Avg. base quality:\t#{@config[:baseq]}"
      puts "Maping quality:\t#{@config[:mapq]}"
      puts "Read group:\twithin #{@config[:within]}bp"
      puts "paired B and F:\twithin #{@config[:pair_within]}bp"
      puts "Support reads for making consensus sequence:\t#{@config[:support_reads]}"
      puts "mimimum clipping fragment base:\t#{@config[:clip_length]}bp"
      puts "support clip length:\t#{@config[:support_clip_length]}bp"
      puts "bam:\t#{@config[:bam]}"
      puts "chr:\t#{@config[:chr]}"
      puts "outd:\t#{@config[:outd]}"
      puts "indelsize:\t#{@config[:indelsize]}"
      puts "reffa:\t#{@config[:reffa]}"
      puts "glsearch:\t#{@config[:glsearch]}"
      puts "glsearch mat:\t#{@config[:glsearch_mat]}"
      puts "mafft:\t#{@config[:mafft]}"
      puts "samtools:\t#{@config[:samtools]}"
      puts "temp:\t#{@config[:temp]}"
      puts "thread:\t#{@config[:thread]}"
      puts "exclude-region:\t#{@config[:exclude_region]}"
      puts
    end

    def default_config
      {
        baseq:               20, # avg base quality in clipping seqs or in unmapped reads
        mapq:                20, # mapping quality
        within:              3,  # readをgroupingする際のbase数: 5bp以内をgrouping
        pair_within:         5,  # LIのpairを作る時: 5bp以内をgrouping
        alt_read_depth:      5,  # alt read depth(non-ref)のminimum数
        support_reads:       3,  # groupingする際のminimum read数
        clip_length:         5,  # clipping fragmentのmimimum base数
        support_clip_length: 5,  # read group clipping fragmentのmimimum base数
        glsearch:            'glsearch36',
        glsearch_mat:        File.join(File.expand_path("", File.dirname(__FILE__)), "..", "data", "mydna.mat"),
        mafft:               'mafft',
        samtools:            'samtools',
        temp:                Dir.tmpdir,
        thread:              1,
      }
    end
  end
end
