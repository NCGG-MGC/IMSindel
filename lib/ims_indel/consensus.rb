require 'bio'
require 'tempfile'
require_relative 'read'
require_relative 'consensus_util'

module IMSIndel
  class Consensus
    include ConsensusUtil

    attr_reader :backward_clip_consensus,
      :forward_clip_consensus,
      :short_indel_consensus,
      :paired_long_insertions,
      :unpaired_long_indels,
      :mafft_inputs

    def initialize(temp_path, thread, mafft, keep_mafft_input)
      @temp_path = temp_path
      @thread = thread
      @mafft = mafft
      @mafft_inputs = {} if keep_mafft_input
    end

    def indel_list
      (short_indel_consensus + paired_long_insertions + unpaired_long_indels).sort do |a,b|
        a.start_pos <=> b.start_pos
      end
    end

    def make_support_read_consensus(read_group, alt_read_depth)
      @backward_clip_consensus = make_consensus_with_threads(read_group.backward_clips, :B, alt_read_depth)
      @forward_clip_consensus  = make_consensus_with_threads(read_group.forward_clips, :F, alt_read_depth)
      @short_indel_consensus   = make_consensus_with_threads(read_group.non_clips, :SID, alt_read_depth)
    end

    def make_consensus_seq(reads, pair_within, alt_read_depth, max_indel_size)
      # >long insertion 
      # find long insertion pairs:
      # [flg_sttpos_endpos_conseq_readcntB, flg_sttpos_endpos_conseq_readcntF]
      insertion_pair = find_insert_pair(@backward_clip_consensus, @forward_clip_consensus, pair_within)

      # detect long insertion pairs:
      # [:ULI, clip_B_stt, clip_B_stt, consensus, ttl_depth]
      # [:LI_wU, clip_stt, clip_end, new_consensus, new_ttl_depth]
      # [:LI, clip_stt, clip_end, consensus, ttl_depth]
      paired_long_insertions = detect_long_insertion(insertion_pair,
                                                     reads.pos2max_backward_clip,
                                                     reads.pos2max_forward_clip,
                                                     reads.unmapped_reads,
                                                     alt_read_depth)

      # >long deletion 
      ins_l_B = insertion_pair.map { |i| i[0] }
      ins_l_F = insertion_pair.map { |i| i[1] }
      del_l_B = @backward_clip_consensus - ins_l_B
      del_l_F = @forward_clip_consensus - ins_l_F
      # update p_indel_l (long indelに加え long deletionを加える)
      # [:LD, clip_chrpos_all.min, clip_chrpos_all.max, consensus, ttl_depth]
      paired_long_insertions, p_del_l_B, p_del_l_F = make_delete_pair(del_l_B, del_l_F, paired_long_insertions, max_indel_size, alt_read_depth)

      # >detection of unpaired indel list
      up_del_l_B = del_l_B - p_del_l_B.uniq
      up_del_l_F = del_l_F - p_del_l_F.uniq
      @unpaired_long_indels = (up_del_l_B + up_del_l_F).find_all{|consensus| consensus.depth >= alt_read_depth} # minimum alt_read_depth check
      @paired_long_insertions = paired_long_insertions
    end


    private

    def make_consensus_with_threads(clips, read_type, alt_read_depth)
      ts = []
      total_size = clips.size
      nthreads = [@thread, total_size].min
      i = 1
      clips.each_slice(clips.size / nthreads) do |slice|
        j = i + slice.size
        ts << Thread.new do
          cons = make_consensus(slice, read_type, alt_read_depth, j, total_size)
          Thread.current['cons'] = cons
        end
        i += slice.size
      end
      ts.each(&:join)
      consensus_all  = []
      ts.each { |t| consensus_all.concat(t['cons']) }
      consensus_all.sort!{|a,b| a.start_pos <=> b.start_pos} # start pos sort
      return consensus_all
    end

    def make_consensus(clips, read_type, alt_read_depth, start_index, total_size)
      consensus_all  = []

      clips.each.with_index(start_index) do |group_reads, index|
        consensus, read_cnt = mafft_consensus(group_reads, 0.8) # %identity = 0.8, 1.0 is too strict
        consensus, trim_flag = trim_consensus(consensus)

        #  consensus = "aaagtgttg?ggtagaaggaaaggaaggaaagagagaaggggaaggaactg?gaga"
        if consensus.empty? || consensus.count("?") > 5 # SNP check and sequence error, if the number of SNPs is more than 5 ( > 5)
          next
        end
        # SNP check, if the number of SNPs is less than 3 ( <= 2)
        consensus, read_cnt = mafft_consensus(group_reads, 0.5) # %identity = 0.5, consensus updated
        if read_type == :SID && read_cnt < alt_read_depth # short indel depth < 5はskip
          next
        end

        consensus = trim_consensus_with_flag(consensus, trim_flag)
        clip_chrpos_min = group_reads[0].start_pos < group_reads[0].end_pos ? group_reads[0].start_pos : group_reads[0].end_pos
        clip_chrpos_max = group_reads[0].start_pos > group_reads[0].end_pos ? group_reads[0].start_pos : group_reads[0].end_pos

        group_reads.each do |read|
          if read.start_pos < read.end_pos
            clip_chrpos_min = read.start_pos if read.start_pos < clip_chrpos_min
            clip_chrpos_max = read.end_pos if read.end_pos > clip_chrpos_max
          else
            clip_chrpos_min = read.end_pos if read.end_pos < clip_chrpos_min
            clip_chrpos_max = read.start_pos if read.start_pos > clip_chrpos_max
          end
        end
        consensus_read = Read.new(type: read_type, start_pos: clip_chrpos_min, end_pos: clip_chrpos_max, seq: consensus, depth: read_cnt)
        consensus_all << consensus_read
        if @mafft_inputs
          @mafft_inputs[consensus_read] = group_reads
        end
      end
      return consensus_all
    end

    # making a consensus seq using MAFFT
    # input1: reads = [Read(map_status, 16161113, 16161113, "CCCAGACGTGTAGAGCTTATGAAAAAAAAAAAAAAAAAAAAAAAA"),
    #                  Read(map_status, 16161114, 16161114, "CCCAGACGTGTAGAGCTTATGAAAAAAAAAAAAAA")]
    # input2: percentID = 0.6
    # input3: outd = tmp/mafft/B_fragment
    def mafft_consensus(reads, percentID)
      tmp = Tempfile.new("maffttmp", @temp_path)
      reads.each.with_index(1) do |read_inf, index|
        tmp.puts ">#{read_inf.type}_#{read_inf.start_pos}_#{read_inf.end_pos}-v#{index}"
        tmp.puts read_inf.seq.upcase
      end
      tmp.close
      td = @temp_path ? "TMPDIR=#{@temp_path}" : ""
      res =`#{td} #{@mafft} --nuc --ep 0.0 --op 1 --genafpair --maxiterate 1000 #{tmp.path} 2>/dev/null`
      tmp.close(true)
      unless $?.success?
        raise "mafft exec error: #{$?}"
      end

      # makeing a consensus seq
      align_reads = {}
      res.split("\n>").each do |align_read|
        align_read_ary = align_read.split("\n")
        if align_read_ary.last == ">"
          if align_read_ary[0].start_with?('>')
            read_name = align_read_ary[0][1..-1]
          else
            read_name = align_read_ary[0]
          end
          align_reads[read_name] = align_read_ary[1..-2].join("")
        else
          read_name = align_read_ary[0]
          align_reads[read_name] = align_read_ary[1..-1].join("")
        end
      end

      aln = Bio::Alignment.new(align_reads.values.sort)
      align_reads_names = []
      consensus = aln.consensus_string(percentID, gap_mode: -1) # threshold =%id

      # tcctcgtgg---tcggctaact-------------------------------------------------------  B_136582615_136582615-v90
      # tcctcgtgg---tcggctaactcctgcaaagcctgagtattctttcatttcatggtgagttttaaatt---------  B_136582615_136582615-v91
      # tcctcgtgg---tcggctaactcctgcaaagcctgagtattctttcatttcatggtgagttttaaatt---------  B_136582615_136582615-v91
      # tcctcgtggAGGtcggctaactcctgcaaagcctgagtattctttcatttcatggtgagttttaaatt---------  B_136582615_136582615-v91
      check = Hash.new(0) # depth1の場所を探し、trimする
      align_reads.each do |read_name, align_seq|
        read_name = read_name[1..-1] if read_name.start_with?(">")
        align_seq.each_char.with_index{ |allele, num| check[num] += 1 if allele != "-" }
        align_reads_names << [align_seq, read_name]
      end
      max_num = check.keys.max

      new_cons = []
      if align_reads_names.size > 2 # multiple-alignmentの場合
        # tcctcgtgg---tcggctaact-------------------------------------------------------  B_136582615_136582615-v90
        # ---tcgtgg---tcggctaactcctgcaaagcctgagtattctttcatttcatggtgagttttaaatt---------  B_136582615_136582615-v91
        # 最初の数文字と最後の数文字はdepth1でも消さない
        # >最初
        bef_index = -1
        flg = 0
        check.sort_by { |k, v| k }.each do |index, cnt|
          if flg == 0 and cnt == 1
            bef_index = index
          elsif flg == 1
            break
          else
            flg = 1
          end
        end
        # >最後
        aft_index = max_num + 1
        flg = 0
        check.sort_by{|k,v|k}.reverse.each do |index, cnt|
          if flg == 0 and cnt == 1
            aft_index = index
          elsif flg == 1
            break
          else
            flg = 1
          end
        end

        # align_reads_namesのチェック
        align_reads_names.each do |align_seq, read_name|
          new_align_seq = ""
          align_seq.each_char.with_index do |seq, num|
            if num <= bef_index || aft_index <= num # 最初と最後のdepth1
              new_align_seq += seq 
            elsif check[num] != 1
              new_align_seq += seq 
            end
          end
        end
        consensus.each_char.with_index do |seq, num|
          if num <= bef_index or aft_index <= num # 最初と最後のdepth1
            new_cons << seq 
          elsif check[num] != 1  
            new_cons << seq 
          end
        end

        # pairwise-alignmentのときは特になにもせずO.K.
      else
        new_cons = [consensus]
      end
      new_cons = new_cons.join("")

      return new_cons, reads.size
    end

    # find insertion pair from cons_Bs and cons_Fs
    def find_insert_pair(backward_consensus, forward_consensus, pair_within)
      insert_pair = []

      forward_consensus_dup = forward_consensus.dup
      backward_consensus.each do |consensus_b|
        start_pos_b = consensus_b.start_pos - pair_within
        end_pos_b = consensus_b.end_pos + pair_within

        skip = 0
        forward_consensus_dup.each_with_index do |consensus_f, index|
          start_pos_f = consensus_f.start_pos - pair_within
          end_pos_f = consensus_f.end_pos + pair_within
          if end_pos_f < start_pos_b
            skip = index
            next
          elsif (start_pos_f <= end_pos_b && end_pos_b <= end_pos_f) || (start_pos_f <= start_pos_b && start_pos_b <= end_pos_f)
            insert_pair << [consensus_b, consensus_f]
          elsif end_pos_b < start_pos_f
            forward_consensus_dup = forward_consensus_dup[skip..-1] # update
            break
          end
        end
      end
      return insert_pair
    end

    # from insert_pairs, find complete insertion and uncomplete insertion
    def detect_long_insertion(insert_pair, clip_Bs2, clip_Fs2, unmap_seq, alt_read_depth=5)
      paired_indel_list = []
      insert_pair.each do |consensus_b, consensus_f|
        total_depth = consensus_b.depth + consensus_f.depth
        clip_chrpos_all = [consensus_b.start_pos, consensus_b.end_pos, consensus_f.start_pos, consensus_f.end_pos]
        clip_start = clip_chrpos_all.min
        clip_end = clip_chrpos_all.max
        consensus, _ = mafft_consensus([consensus_b, consensus_f], 1.0) # %identity = 1.0
        consensus, trim_flag = trim_consensus(consensus)
        # ---------------------------------------------------------------------------

        if consensus.count("?") > 2 # SNP check, if the number of SNPs is more than 2 ( > 2), uncomplete Long Insertion
          #puts "uncomplete long insertion..."

          # re-alignement using unmapped read
          can_use_unmaps = unmap_seq.find_all do | unmap|
            (unmap.start_pos < clip_start && clip_start < unmap.end_pos) || (unmap.start_pos < clip_end && clip_end < unmap.end_pos)
          end
          if can_use_unmaps.empty? # unmap read is nothing
            clip_B_stt = consensus_b.start_pos # clip_B start pos
            clip_F_end = consensus_f.end_pos # clip_F end pos
            upper_LI  = clip_Bs2[clip_B_stt]
            bottom_LI = clip_Fs2[clip_F_end]
            consensus = "#{upper_LI}-----#{bottom_LI}"
            if total_depth >= alt_read_depth # end pos is equal to sttpos
              paired_indel_list << Read.new(type: :ULI, start_pos: clip_B_stt, end_pos: clip_B_stt, seq: consensus, depth: total_depth)
              @mafft_inputs[paired_indel_list.last] = [consensus_b, consensus_f] if @mafft_inputs
            end

          else # unmap reads exist
            new_group_reads = ([consensus_b, consensus_f] + can_use_unmaps)
            new_consensus, _ = mafft_consensus(new_group_reads, 1.0) # make a consensus seq with the unmapped reads
            new_consensus, trim_flag = trim_consensus(new_consensus)

            if new_consensus.empty? || new_consensus.count("?") > 2 # SNP check, if the number of SNPs is more than 2 ( >=3)
              clip_B_stt = consensus_b.start_pos # clip_B start pos
              clip_F_end = consensus_f.end_pos # clip_F end pos
              upper_LI  = clip_Bs2[clip_B_stt]
              bottom_LI = clip_Fs2[clip_F_end]
              consensus = "#{upper_LI}-----#{bottom_LI}"
              if total_depth >= alt_read_depth # end pos is equal to sttpos
                paired_indel_list << Read.new(type: :ULI, start_pos: clip_B_stt, end_pos: clip_B_stt, seq: consensus, depth: total_depth)
                @mafft_inputs[paired_indel_list.last] = new_group_reads if @mafft_inputs
              end
            else 
              new_consensus, _ = mafft_consensus(new_group_reads, 0.5) # make a consensus seq with the unmapped reads
              new_consensus = trim_consensus_with_flag(new_consensus, trim_flag)

              new_total_depth = new_group_reads.inject(0){|res, read| res += read.depth}
              if new_total_depth >= alt_read_depth ####
                paired_indel_list << Read.new(type: :LI_wU, start_pos: clip_start, end_pos: clip_end, seq: new_consensus, depth: new_total_depth)
                @mafft_inputs[paired_indel_list.last] = new_group_reads if @mafft_inputs
              end
            end
          end

        else
          consensus, _ = mafft_consensus([consensus_b, consensus_f], 0.5) # %identity = 0.5, consensus update
          consensus = trim_consensus_with_flag(consensus, trim_flag)

          #puts "complete long insertion..."
          #puts align_reads_names
          #puts
          if total_depth >= alt_read_depth ####
            paired_indel_list << Read.new(type: :LI, start_pos: clip_start, end_pos: clip_end, seq: consensus, depth: total_depth)
            @mafft_inputs[paired_indel_list.last] = [consensus_b, consensus_f] if @mafft_inputs
          end
        end
      end
      return paired_indel_list
    end

    def make_delete_pair(del_list_B, del_list_F, paired_indel_list, i_size, alt_read_depth=5)
      # find deletion pairs
      del_list_F_dup = del_list_F.dup
      paired_del_list_B = []
      paired_del_list_F = []
      if del_list_B.empty? || del_list_F.empty?
        return paired_indel_list, paired_del_list_B, paired_del_list_F
      end

      del_list_B.each do |readB|
        skip = 0
        del_list_F_dup.each_with_index do |readF, index|
          if readF.end_pos < readB.start_pos
            skip = index
            next
          end
          break if (readB.end_pos + i_size) < readF.start_pos # indel size == i_size:5000

          group_reads = [readB, readF]
          consensus, _ = mafft_consensus(group_reads, 1.0)
          consensus, trim_flag = trim_consensus(consensus)

          if consensus.empty? || consensus.count("?") > 2
            next
          else
            consensus, _ = mafft_consensus(group_reads, 0.5) # update, consensus
            consensus = trim_consensus_with_flag(consensus, trim_flag)

            clip_chrpos_all = [readB.start_pos, readB.end_pos, readF.start_pos, readF.end_pos]
            ttl_depth = readB.depth + readF.depth
            if ttl_depth >= alt_read_depth
              paired_indel_list << Read.new(type: :LD, start_pos: clip_chrpos_all.min,
                                            end_pos: clip_chrpos_all.max, seq: consensus, depth: ttl_depth)

              @mafft_inputs[paired_indel_list.last] = group_reads if @mafft_inputs
            end
            paired_del_list_B << readB
            paired_del_list_F << readF
            #puts "paired long deletion..."
            del_list_F_dup = del_list_F_dup[skip..-1] # update
            skip = 0
          end
        end
      end
      return paired_indel_list, paired_del_list_B, paired_del_list_F
    end
  end
end
