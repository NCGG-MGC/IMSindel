require_relative 'pairwise_consensus'

module IMSIndel
  class IndelDetector

    def initialize(glsearch_bin, glsearch_mat, samtools, temp_path)
      @glsearch_bin = glsearch_bin
      @glsearch_mat = glsearch_mat
      @samtools = samtools
      @temp_path = temp_path
    end

    def indel_call(indel_list, target_chr, max_indel_size, reffa, bam_path, mapq, alt_read_depth, outf)
      all_pos, fin_indels = detect_indels(indel_list, target_chr, max_indel_size, reffa)

      puts "#all_pos: #{all_pos.size}"
      puts "#fin_indels: #{fin_indels.size}"
      fin_indels.default_proc = nil

      if all_pos.empty?
        puts "empty all_pos"
        return
      end
      all_pos_depth = count_depth(all_pos, bam_path, target_chr, mapq)
      puts "#all_pos_depth #{all_pos_depth.size}"

      final_indel_call(fin_indels, all_pos_depth, outf, target_chr, alt_read_depth)
    end

    private

    # indelの開始終了情報 struct
    IndelInfo = Struct.new(:start, :end, :str, :is_same, :key)
    # 見つかったindel struct
    IndelResult = Struct.new(:indel_type, :flag, :depth, :start_pos, :seq)

    # compared with reference seq, detect indels precisely
    def detect_indels(indel_list, target_chr, max_indel_size, reffa)
      @fin_indels = Hash.new {|h,k| h[k] = [] } # {31973646_31973663_tgttatttataataatga" => ["INS_SID_18(depth)_indelsttpos_seq", ..
      @is_checked_pos = {} # same Indel check (temprary use)
      @all_pos    = []
      pairwise_consensus = PairwiseConsensus.new(@glsearch_bin, @glsearch_mat, @temp_path)

      ttl_indel = indel_list.size
      indel_list.each_with_index do |read, index| # read = [":LI(or :ULI or :LD or :B or :F or :SID), clip_chrpos_all.min, clip_chrpos_all.max, consensus, ttl_depth]
        flag = read.type
        sttpos = read.start_pos
        endpos = read.end_pos
        seq = read.seq
        depth = read.depth

        if flag == :ULI # uncomplete long indel
          key = "#{sttpos}_#{endpos}"
          if @is_checked_pos[key]
            #puts "#{index+1}/#{ttl_indel}\t#{flag}_#{sttpos}-#{endpos}\t#{seq}\t#{seq.size-5}+\tdepth: #{depth}: INS\tSAME\n"
            @fin_indels[key] << IndelResult.new("INS", flag, depth, sttpos, seq)
            @all_pos << sttpos
          else
            @is_checked_pos[key] = true
            #puts "#{index+1}/#{ttl_indel}\t#{flag}_#{sttpos}-#{endpos}\t#{seq}\t#{seq.size-5}+\tdepth: #{depth}: INS\n"
            @fin_indels[key] << IndelResult.new("INS", flag, depth, sttpos, seq)
            @all_pos << sttpos
          end
          next
        end

        # cut a part of reference seq
        if sttpos - (max_indel_size/2) < 1
          ref_stt = 1
        else
          ref_stt = sttpos - (max_indel_size/2)
        end
        ref_end = endpos + (max_indel_size/2)
        ref_seq = `samtools faidx #{reffa} #{target_chr}:#{ref_stt}-#{ref_end}`.chomp.split("\n")[1..-1].join("").upcase
        ref  = [:R, ref_stt, ref_end, ref_seq]
        ipt_read = [read.type, read.start_pos, read.end_pos, read.seq]
        group_reads = [ipt_read, ref]
        indel_stt = ipt_read[1] # indel start position before comparison of reference genome, #######################
        consensus, align_reads_names, trim_flag = pairwise_consensus.make_consensus(group_reads, 1.0, nil)
        next if consensus.nil?

        if consensus.empty? || consensus.count("?") > 2 # ?2個はO.K.
          #align_reads_names.each_with_index do |line, subindex|
            #if subindex == 2
              #puts "#{index+1}/#{ttl_indel}\t#{line}\tNOT!!"
            #elsif /^R_/ =~ line.chomp.split("\t")[-1]
              #puts "Ref\t#{line}\n"
            #else
              #puts "Seq\t#{line}\n"
            #end
          #end
          next
        end

        consensus, align_reads_names = pairwise_consensus.make_consensus(group_reads, 0.5, trim_flag)
        next if consensus.nil?

        read_seq, _ = align_reads_names[0].split("\t")
        ref_seq, ref_name   = align_reads_names[1].split("\t")

        # special for LI, 150525
        # ref_seq = cagaaaaa--aaaaa---------ccagtcccgttctaaccagtcccgttctaa, :first --はmiss 
        #puts "ref\n#{ref_seq}\nread\n#{read_seq}\n"

        if (flag == :LI || flag == :LI_wU) && /(^\w+)(-+)(\w+)(-+)(\w+)$/ =~ ref_seq
          #puts "LI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n#{flag}\n"
          fir_str = $1.size
          fir_gap = $2.size
          sec_str = $3.size
          sec_gap = $4.size
          next unless /(^\-+)(\w+)(\-+)$/ =~ read_seq
          indel_info = process_long_insersion(fir_gap, fir_str, sec_gap, sec_str, read_seq, ref_name, ref_seq)
          @fin_indels[indel_info.key] << IndelResult.new("INS", flag, depth, indel_stt, indel_info.str)
          @all_pos << indel_stt
          #puts_indel_info("INS", ref_seq, read_seq, index, ttl_indel, align_reads_names, flag, depth, indel_info)

        elsif /(^\w+)(-+)(\w+)$/ =~ ref_seq # ref = cagaaaaa-----ccagtcccgttctaaccagtcccgttctaa
          # insertion: 右寄せのパターン or 真ん中が
          # 1. mono polymer
          # read = "TTTaaaaGGGGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
          # ref  = "---aaaaGGG--aaaaaaaaaaaaaaaaaaaaaa-------"
          # 2. multi polymer (右寄せ)
          # ref  = "TTTaaaa<a><b> -  - TTTT"
          # read = "---aaaa<a><b><a><b>TTTT"
          fir_str = $1.size
          fir_gap = $2.size
          sec_str = $3.size
          next unless /(^\-+)(\w+)(\-+)$/ =~ read_seq

          if no_overlap?(ref_seq, read_seq, fir_str, fir_gap)
            #puts "Ref\t#{ref_seq}\nSeq\t#{read_seq}"
            #puts "#{index+1}/#{ttl_indel}\tNO Overlap(Ins)!!\n"
            next
          end
          if few_over_lap_insertion?(read_seq, fir_str, fir_gap)
            #puts "Ref\t#{ref_seq}\nSeq\t#{read_seq}"
            #puts "#{index+1}/#{ttl_indel}\tFew Overlap(Ins)!!"
            next
          end
          indel_info = process_right_insersion(fir_gap, fir_str, sec_gap, sec_str, read_seq, ref_name, ref_seq)
          @fin_indels[indel_info.key] << IndelResult.new("INS",  flag, depth, indel_stt, indel_info.str)
          @all_pos << indel_stt
          #puts_indel_info("INS", ref_seq, read_seq, index, ttl_indel, align_reads_names, flag, depth, indel_info)

        elsif ref_seq.count("-") == 0
          # deletion: 右寄せのパターンが多い
          # 1. mono polymer
          # ref  = "TTTaaaaGGGGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
          # read = "---aaaaGGG--aaaaaaaaaaaaaaaaaaaaaa-------"
          # 2. multi polymer
          # ref  = "TTTaaaaATGAGTAGTaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
          # read = "---aaaaAGTAGT---aaaaaaaaaaaaaaaaaaaaaa-------"
          #
          # first str >= 10 and sec str >= 10
          if /(^\-+)(\w+)(\-+)(\w+)(?:\-+)$/ =~ read_seq
            fir_gap  = $1.size
            fir_str  = $2.size
            sec_gap  = $3.size
            sec_str  = $4.size
            if few_over_lap_deletion?(fir_str, sec_str)
              #puts "Ref\t#{ref_seq}\nSeq\t#{read_seq}"
              #puts "#{index+1}/#{ttl_indel}\tFew Overlap(Del)!!\n"
              next
            end
            indel_info = process_deletion(fir_gap, fir_str, sec_gap, sec_str, read_seq, ref_name, ref_seq)
            @fin_indels[indel_info.key] << IndelResult.new("DEL", flag, depth, indel_stt, indel_info.str)
            @all_pos << indel_stt
          else
            #puts "Ref\t#{ref_seq}\nSeq\t#{read_seq}"
            #puts "#{index+1}/#{ttl_indel}\t#{align_reads_names[-1]}"
            #puts "#{index+1}/#{ttl_indel}\t#{flag}\tNOT_INDEL\n"
          end
        end
      end
      @is_checked_pos = nil
      @all_pos.sort!
      @all_pos.uniq!
      return @all_pos, @fin_indels
    end

    # count refeerence allele depth
    # count total depth
    # ------------------------------------------
    def count_depth(all_pos, bam_path, target_chr, mapq)
      pos_max = all_pos.max 
      pos_min = all_pos.min
      all_pos_depth = Hash.new(0)

      SAMReader.new(@samtools, bam_path).each(chr: target_chr,
                                         filter_flag: SAMReader::DUPLICATE,
                                         output_flag: SAMReader::PROPER_PAIR) do |line|
        read_seq = line.read_seq
        next if target_chr != line.chr 
        next if line.map_score <= mapq # Condition1
        next if read_seq.count("N") > 0
        map_status = line.map_status

        # complete match
        if /^(\d+)M$/ =~ map_status # 161M
          match_size = $1.to_i  #  obj_SAM.map_status --> 70M or 101M or

          # B fragment (match + clip seq)
        elsif /^(\d+)M(?:\d+)S$/ =~ map_status # 80M20S
          len_M = $1.to_i
          match_size = len_M
          # B fragment (short INS + clip seq)
        elsif /^(\d+)M(?:\d+)I(\d+)M(?:\d+)S$/ =~ map_status # 70M2I8M20S
          len_M1 = $1.to_i #70
          len_M2 = $2.to_i #8
          match_size = len_M1 + len_M2
          # B fragment (short DEL + clip seq)
        elsif /^(\d+)M(\d+)D(\d+)M(?:\d+)S$/ =~ map_status # 70M2D10M20S
          len_M1 = $1.to_i #70
          len_D  = $2.to_i #2
          len_M2 = $3.to_i #8
          match_size = len_M1 + len_D + len_M2

          # F fragment (clip seq + match)
        elsif /^(?:\d+)S(\d+)M$/ =~ map_status # 20S80M
          len_M = $1.to_i
          match_size = len_M
          # F fragment (clip seq + short INS)
        elsif /^(?:\d+)S(\d+)M(?:\d+)I(\d+)M$/ =~ map_status # 20S70M2I8M
          len_M1 = $1.to_i #70
          len_M2 = $2.to_i #8
          match_size = len_M1 + len_M2
          # F fragment (clip seq + short DEL)
        elsif /^(?:\d+)S(\d+)M(\d+)D(\d+)M$/ =~ map_status # 20S70M2D10M
          len_M1 = $1.to_i #70
          len_D  = $2.to_i #2
          len_M2 = $3.to_i #8
          match_size = len_M1 + len_D + len_M2

          # -----------------------------------------
          # short INDEL
          # -----------------------------------------
          # 1 insertion 
        elsif /^(\d+)M(\d+)I(\d+)M$/ =~ map_status # 70M2I28M
          len_M1 = $1.to_i #70
          #len_I  = $2.to_i #2
          len_M2 = $3.to_i #28
          match_size = len_M1 + len_M2
          # 1 deletion
        elsif /^(\d+)M(\d+)D(\d+)M$/ =~ map_status # 70M5D30M
          len_M1 = $1.to_i #70
          len_D  = $2.to_i #5
          len_M2 = $3.to_i #30
          match_size = len_M1 + len_D + len_M2

          # 1 del + 1 del
        elsif /^(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M$/ =~ map_status # 70M 2D 10M 5D 28M
          len_M1 = $1.to_i #70
          len_D1 = $2.to_i #2
          len_M2 = $3.to_i #10
          len_D2 = $4.to_i #5
          len_M3 = $5.to_i #28
          match_size  = len_M1 + len_D1 + len_M2 + len_D2 + len_M3
          # 1 del + 1 ins
        elsif /^(\d+)M(\d+)D(\d+)M(?:\d+)I(\d+)M$/ =~ map_status # 70M 2D 10M 5I 28M
          len_M1 = $1.to_i #70
          len_D1 = $2.to_i #2
          len_M2 = $3.to_i #10
          len_M3 = $4.to_i #28
          match_size  = len_M1 + len_D1 + len_M2 + len_M3
          # 1 ins + 1 del
        elsif /^(\d+)M(?:\d+)I(\d+)M(\d+)D(\d+)M$/ =~ map_status # 70M 2I 10M 5D 28M
          len_M1 = $1.to_i #70
          len_M2 = $2.to_i #10
          len_D1 = $3.to_i #5
          len_M3 = $4.to_i #28
          match_size  = len_M1 + len_M2 + len_D1 + len_M3
          # 1 ins + 1 ins
        elsif /^(\d+)M(?:\d+)I(\d+)M(?:\d+)I(\d+)M$/ =~ map_status # 70M 2I 10M 5I 28M
          len_M1 = $1.to_i #70
          len_M2 = $2.to_i #10
          len_M3 = $3.to_i #28
          match_size  = len_M1 + len_M2 + len_M3
        else
          next
        end

        sttpos = line.chrpos
        endpos = sttpos + match_size - 1
        next if endpos < pos_min
        break if pos_max < sttpos

        skip = 0
        all_pos.each_with_index do |pos, index|
          if pos < sttpos
            skip = index
            next
          elsif sttpos <= pos && pos <= endpos
            all_pos_depth[pos] += 1  
          elsif endpos < pos
            all_pos = all_pos[skip..-1] # all_pos update
            break 
          end
        end
      end
      return all_pos_depth
    end

    def final_indel_call(fin_indels, all_pos_depth, outf, iptchr, alt_read_depth=5)
      f = open(outf, "w")
      f.puts ">indel_type\tcall_type\tchr\tsttpos\tendpos\tindel_length\tindel_str\t#indel_depth\t#ttl_depth\tdetails(indelcall_indeltype_depth)\tclip_sttpos\tdepth(>=10)"
      
      fin_indels.sort_by{|k,v| k.split("_")[0].to_i}.each do |sttpos_endpos, indel_results| # start pos sort,  deppth_ary_str =["INS_SID_18(depth)_indelsttpos_seq", ..
        sttpos, endpos = sttpos_endpos.split("_")
        indel_call  = ''
        #
        # 複数種類のindelがある場合がある(それぞれalt depth>=5), +AAA/+AA
        # 
        indelstr_check = []
        indel_results.each do |indel_result| 
          indelstr_check << indel_result.seq if indel_result.depth >= alt_read_depth #next alt allele read count < 5
        end
        indelstr_check.uniq!
        next if indelstr_check.empty?
        
        tmp_indellen = []
        tmp_indelstr = []
        tmp_depthary = []
        tmp_clippos  = []
        tmp_depth    = []
        tmp_ttlcnt   = []
        tmp_freq     = []
        tmp_calls    = []
        indelstr_check.each do |fcs_indel|
          depth_ary    = []
          dep_clippos  = {}
          dep_indelstr = {}
          type_depth   = Hash.new{|h,k| h[k] = []}
          ins_del      = ''
          
          indel_results.each do |indel_result|
            #IndelResult(:indel_type, :flag, :depth, :start_pos, :seq)
            if indel_result.seq == fcs_indel
              ins_del = indel_result.indel_type  #indel_call
              depth_ary << "#{indel_result.indel_type}_#{indel_result.flag}_#{indel_result.depth}_#{indel_result.start_pos}"
              dep_clippos[indel_result.depth]  = indel_result.start_pos # {dep:36 => clippos:1927, 5=>1927, 
              dep_indelstr[indel_result.depth] = indel_result.seq       # {dep:36 => AAAAA, 
              type_depth[indel_result.flag] << indel_result.depth       #  {SID => [dep:36, dep:10], F =>5
            end
          end
          
          # indel type別のmax depthを足してtotal depthのカウント
          depth = type_depth.values.inject(0){|sum, each_depths| sum += each_depths.max }

          # total depth
          ttlcnt_clippos_indelstr = {}
          all_depth = type_depth.values.flatten # [10, 20, 30]
          all_depth.each do |each_depth|
            each_clippos   = dep_clippos[each_depth]
            each_indel_str = dep_indelstr[each_depth]
            ttl_cnt        = all_pos_depth[each_clippos]
            ttlcnt_clippos_indelstr[ttl_cnt] = [each_clippos, each_indel_str]
          end
          ttl_cnt = ttlcnt_clippos_indelstr.keys.max # maximam
          clippos, indel_str = ttlcnt_clippos_indelstr[ttl_cnt]
          ttl_cnt = depth if depth > ttl_cnt # fix
          freq    = depth / ttl_cnt.to_f
          
          # 保存
          tmp_indellen << indel_str.size
          tmp_indelstr << indel_str
          tmp_depth    << depth
          tmp_ttlcnt   << ttl_cnt
          tmp_depthary << depth_ary.join(" ")
          tmp_clippos  << clippos
          tmp_freq     << freq
          tmp_calls    << ins_del
        end
        

        if tmp_indelstr.size == 1
          indelsize  = tmp_indellen.first
          indel_str  = tmp_indelstr.first
          depth      = tmp_depth.first
          ttl_cnt    = tmp_ttlcnt.first
          depth_ary  = tmp_depthary.first
          clippos    = tmp_clippos.first
          freq       = tmp_freq.first
          indel_call = tmp_calls.first

          if 0.70 <= freq # homo
            if ttl_cnt >= 10
              f.puts "#{indel_call}\tHomo\t#{iptchr}\t#{sttpos}\t#{endpos}\t#{indelsize}\t#{indel_str}\t#{depth}\t#{ttl_cnt}\t#{depth_ary}\t#{clippos}\tHigh"
            else
              f.puts "#{indel_call}\tHomo\t#{iptchr}\t#{sttpos}\t#{endpos}\t#{indelsize}\t#{indel_str}\t#{depth}\t#{ttl_cnt}\t#{depth_ary}\t#{clippos}\tLow"
            end

          elsif (0.15 < freq && freq <= 0.70) # hetero
            if ttl_cnt >= 10
              f.puts "#{indel_call}\tHete\t#{iptchr}\t#{sttpos}\t#{endpos}\t#{indelsize}\t#{indel_str}\t#{depth}\t#{ttl_cnt}\t#{depth_ary}\t#{clippos}\tHigh"
            else
              f.puts "#{indel_call}\tHete\t#{iptchr}\t#{sttpos}\t#{endpos}\t#{indelsize}\t#{indel_str}\t#{depth}\t#{ttl_cnt}\t#{depth_ary}\t#{clippos}\tLow"
            end

          else # low frequency
            if ttl_cnt >= 10
              f.puts "#{indel_call}\tHete\t#{iptchr}\t#{sttpos}\t#{endpos}\t#{indelsize}\t#{indel_str}\t#{depth}\t#{ttl_cnt}\t#{depth_ary}\t#{clippos}\tHigh,Lowfreq"
            else
              f.puts "#{indel_call}\tHete\t#{iptchr}\t#{sttpos}\t#{endpos}\t#{indelsize}\t#{indel_str}\t#{depth}\t#{ttl_cnt}\t#{depth_ary}\t#{clippos}\tLow,Lowfreq"
            end
          end


        #  それぞれ違うもののhetero. c.g. +AA/+AAA
        elsif tmp_indelstr.size == 2
          cnt = 0
          tmp_ttlcnt.each do |ttl|
            cnt += 1 if ttl.to_i >= 10
          end
          if cnt == 2
            f.puts "#{tmp_calls.join(",")}\tHete\t#{iptchr}\t#{sttpos}\t#{endpos}\t#{tmp_indellen.join(",")}\t#{tmp_indelstr.join(",")}\t#{tmp_depth.join(",")}\t#{tmp_ttlcnt.join(",")}\t#{tmp_depthary.join(",")}\t#{tmp_clippos.join(",")}\tHigh" 
          else
            f.puts "#{tmp_calls.join(",")}\tHete\t#{iptchr}\t#{sttpos}\t#{endpos}\t#{tmp_indellen.join(",")}\t#{tmp_indelstr.join(",")}\t#{tmp_depth.join(",")}\t#{tmp_ttlcnt.join(",")}\t#{tmp_depthary.join(",")}\t#{tmp_clippos.join(",")}\tLow" 
          end
        end
      end
      f.close
    end

    # indel情報print
    # 行数取るのでメソッド化している
    def puts_indel_info(type, ref_seq, read_seq, index, ttl_indel, align_reads_names, flag, depth, indel_info)
      puts "Ref\t#{ref_seq}"
      puts "Seq\t#{read_seq}"
      puts "#{index+1}/#{ttl_indel}\t#{align_reads_names[-1]}"
      print "#{index+1}/#{ttl_indel}\t#{flag}_#{indel_info.start}-#{indel_info.end}\t#{indel_info.str}\t#{indel_info.str.size}\tdepth: #{depth}: #{type}"
      print "\tSAME" if indel_info.is_same
      puts "\n"
    end

    # long insertionのstart, end posを計算
    # (flag == :LI || flag == :LI_wU) && /(^\w+)(-+)(\w+)(-+)(\w+)$/ =~ ref_seqのものが対象
    def process_long_insersion(fir_gap, fir_str, sec_gap, sec_str, read_seq, ref_name, ref_seq)
      is_same = true
      ref_stt = ref_name.split("_")[1].to_i

      if fir_gap > sec_gap # fir_gap is LI
        ins_str = read_seq[fir_str..(fir_str + fir_gap - 1)].upcase
        ins_stt = ref_stt + fir_str - 1
      else
        ins_str = read_seq[(fir_str + fir_gap + sec_str)..(fir_str + fir_gap + sec_str + sec_gap - 1)].upcase
        ins_stt = ref_stt + fir_str + sec_str - 1
      end
      ins_end = ins_stt
      key = "#{ins_stt}_#{ins_end}"
      unless @is_checked_pos[key]
        @is_checked_pos[key] = true
        is_same = false
      end
      return IndelInfo.new(ins_stt, ins_end, ins_str, is_same, key)
    end

    # insertionのオーバーラップ判定
    def no_overlap?(ref_seq, read_seq, fir_str, fir_gap)
      # ---------------------
      # NO indel
      # ---------------------
      # ref AAAA-----GGGGGG
      # seq ----AGAGAGGGG--
      if read_seq[0..(fir_str-1)].split("").uniq.size == 1
        return true
      end
      # ref AAA---GGGGGG
      # seq -AAAGA------
      if read_seq[(fir_str + fir_gap)..-1].split("").uniq.size == 1
        return true
      end
      return false
    end

    # insertionのオーバーラップ判定。すくなすぎるとだめ
    def few_over_lap_insertion?(read_seq, fir_str, fir_gap)
      # ------------------
      # overlap < 5はダメ(insertion)
      # ------------------               
      # ref AAAAA---GGGGGG
      # seq ---AAAAAGGG---

      # before overlap
      before     = read_seq[0..(fir_str-1)].size 
      before_gap = read_seq[0..(fir_str-1)].count("-")
      # after overlap
      after     = read_seq[(fir_str + fir_gap)..-1].size
      after_gap = read_seq[(fir_str + fir_gap)..-1].count("-")
      if before - before_gap < 5 || after - after_gap < 5
        return true
      end
      return false
    end

    # deletionのオーバーラップ判定。少なすぎるとだめ
    def few_over_lap_deletion?(fir_str, sec_str)
      # overlap < 10はダメ(deletion)
      # ref : AAAAAGGGGAAAAA
      # seq : --AA--------A-
      return fir_str < 10 || sec_str < 10
    end

    # insertion: 右寄せのパターン or 真ん中が
    # 1. mono polymer
    # read = "TTTaaaaGGGGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    # ref  = "---aaaaGGG--aaaaaaaaaaaaaaaaaaaaaa-------"
    # 2. multi polymer (右寄せ)
    # ref  = "TTTaaaa<a><b> -  - TTTT"
    # read = "---aaaa<a><b><a><b>TTTT"
    def process_right_insersion(fir_gap, fir_str, sec_gap, sec_str, read_seq, ref_name, ref_seq)
      ref_stt = ref_name.split("_")[1].to_i
      ref_str = ref_seq[0..(fir_str - 1)].upcase  # CGGAAA, #######
      ins_str = read_seq[fir_str..(fir_str + fir_gap - 1)].upcase  # AA
      ins_stt = 0
      ins_end = 0
      is_same = true

      # left aligned:
      # inst_str == mono polymer: GGG
      if ins_str.split("").uniq.size == 1
        mono = ins_str[0..0] 
        cnt = 0
        ref_str.split("").reverse.each do |nuc|
          if nuc == mono
            cnt += 1
          else
            break
          end
        end

        ins_stt = ref_stt + fir_str - 1 - cnt
        ins_end = ins_stt

        # left aligned:
        #
        # multi polymer
      else
        ins_size = ins_str.size
        ref_size = ref_str.size
        # ref_seq  = "TTTaaaaATAATA---GTTT"
        # read_seq = "---aaaaATAATAATAGTTTT"
        if ref_str[-ins_size..-1] == ins_str # まず繰り返しがある
          cnt = 0
          iterate = (ref_size / ins_size).to_i
          for i in 1..iterate
            sttpos = ins_size * -i
            endpos = (ins_size * -i) + ins_size - 1
            nuc = ref_str[sttpos..endpos]
            if nuc == ins_str
              cnt += 1
            else
              break
            end
          end

          ins_stt = ref_stt + fir_str - 1 - (cnt * ins_size)
          ins_end = ins_stt
          ref_str = ref_seq[0..((fir_str - 1) - (cnt * ins_size))].upcase  # CGGAAA,
          # ins_strを前半と後半に分けて行う
          for i in 1..ins_size
            aft_str = ins_str[-i..-1]
            bef_str = ins_str[0..-(i+1)]
            if aft_str != ref_str[-i..-1] 
              if i == 1
                break
              elsif i > 1
                aft_str = ins_str[-(i-1)..-1]
                bef_str = ins_str[0..-i]
                ins_str = aft_str + bef_str
                ins_stt = ins_stt - (i-1)
                ins_end = ins_stt
                break
              end
            end
          end

          # 繰り返しがない
        else
          ins_stt = ref_stt + fir_str - 1
          ins_end = ins_stt
          # ins_strを前半と後半に分けて行う
          for i in 1..ins_size
            aft_str = ins_str[-i..-1]
            bef_str = ins_str[0..-(i+1)]

            if aft_str != ref_str[-i..-1] 
              if i == 1
                break
              elsif i > 1
                aft_str = ins_str[-(i-1)..-1]
                bef_str = ins_str[0..-i]
                ins_str = aft_str + bef_str
                ins_stt = ins_stt - (i-1)
                ins_end = ins_stt
                break
              end
            end

          end
        end
      end
      key = "#{ins_stt}_#{ins_end}"
      unless @is_checked_pos[key]
        @is_checked_pos[key] = true
        is_same = false
      end
      return IndelInfo.new(ins_stt, ins_end, ins_str, is_same, key)
    end

    # deletion: 右寄せのパターンが多い
    # 1. mono polymer
    # ref  = "TTTaaaaGGGGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    # read = "---aaaaGGG--aaaaaaaaaaaaaaaaaaaaaa-------"
    # 2. multi polymer
    # ref  = "TTTaaaaATGAGTAGTaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    # read = "---aaaaAGTAGT---aaaaaaaaaaaaaaaaaaaaaa-------"
    #
    # first str >= 10 and sec str >= 10
    def process_deletion(fir_gap, fir_str, sec_gap, sec_str, read_seq, ref_name, ref_seq)
      ref_str = ref_seq[0..(fir_gap + fir_str-1)].upcase #TTTaaaaGGG, or TTTaaaaATGAGT ###############

      del_str = ref_seq[(fir_gap + fir_str)..(fir_gap + fir_str + sec_gap - 1)].upcase 
      ref_stt = ref_name.split("_")[1].to_i
      # left aligned
      # del_str == mono polymer: GGG
      if del_str.split("").uniq.size == 1
        mono = del_str[0..0] #G,
        cnt = 0
        ref_str.split("").reverse.each do |nuc|
          if nuc == mono
            cnt += 1
          else
            break
          end
        end
        del_stt = ref_stt + fir_gap + fir_str - 1 + 1 - cnt
        del_end = ref_stt + fir_gap + fir_str + sec_gap - 1 - cnt

        # del_str == AGT, fragmentの繰り返し
        # ref_str = TTTaaaaATGAGT
      else
        cnt = 0 
        ref_size = ref_str.size
        del_size = del_str.size
        iterate = (ref_size / del_size).to_i
        for i in 1..iterate
          sttpos = del_size * -i
          endpos = (del_size * -i) + del_size - 1
          nuc = ref_str[sttpos..endpos]
          if nuc == del_str
            cnt += 1
          else
            break
          end
        end
        # ref  = "TTTaaaaATGAGTAGTaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
        # read = "---aaaaAGTAGT---aaaaaaaaaaaaaaaaaaaaaa-------"
        del_stt = ref_stt + fir_gap + fir_str - 1 + 1 - (cnt * del_size)
        del_end = ref_stt + fir_gap + fir_str + sec_gap - 1 - (cnt * del_size)
        ref_str = ref_seq[0..(fir_gap + fir_str - 1 - (cnt * del_size))].upcase
        # del_strを前半と後半に分けて行う
        for i in 1..del_size
          aft_str = del_str[-i..-1]
          bef_str = del_str[0..-(i+1)]
          if aft_str != ref_str[-i..-1] 
            if i == 1
              break
            elsif i > 1
              aft_str = del_str[-(i-1)..-1]
              bef_str = del_str[0..-i]
              del_str = aft_str + bef_str

              del_stt = del_stt - (i-1)
              del_end = del_stt + del_str.size - 1
              break
            end
          end
        end
      end
      key = "#{del_stt}_#{del_end}"
      is_same = true
      unless @is_checked_pos[key]
        @is_checked_pos[key] = true
        is_same = false
      end
      return IndelInfo.new(del_stt, del_end, del_str, is_same, key)
    end
  end
end
