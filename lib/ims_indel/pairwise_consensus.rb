require_relative 'consensus_util'

module IMSIndel
  class PairwiseConsensus
    include ConsensusUtil

    def initialize(glsearch_bin, glsearch_mat, temp_path)
      @glsearch_bin = glsearch_bin
      @glsearch_mat = glsearch_mat
      @temp_path = temp_path
    end

    # --------------------------------------------------------------------------------------------------------------------------------------------
    #
    # making a consensus seq using glsearch (pairwise-alignment)
    #
    # input1: reads = [[flg, 16161113, 16161113, "CCCAGACGTGTAGAGCTTATGAAAAAAAAAAAAAAAAAAAAAAAA"], [flg, 16161114, 16161114, "CCCAGACGTGTAGAGCTTATGAAAAAAAAAAAAAA"]]
    # input2: percentID = 0.5 or 1.0
    #
    # --------------------------------------------------------------------------------------------------------------------------------------------
    #
    def make_consensus(reads, percent_id, trim_flag)
      glsearch_result = exec_glsearch(reads)

      gap_priority_alignment     = glsearch_result[:gap_priority_alignment]
      non_gap_priority_alignment = glsearch_result[:non_gap_priority_alignment]
      gap_flags                  = glsearch_result[:gap_flags]
      non_gap_flags              = glsearch_result[:non_gap_flags]
      flag_name                  = glsearch_result[:flag_name]

      # -----------------------------
      # check for gap_priority (gap優先)
      # -----------------------------
      sttpos  = 0
      gap_count = 0
      count = 0
      gap_priority_alignment.split("\n").each do |line|
        count += 1 if line.start_with?('>') # /^>/ =~ line
        case count
        when 0
          next
        when 1
          sttpos = $1.to_i if /^global\/local score.+:(.+)-\d+\)$/ =~ line
          id = line.split("\s").first

          unless gap_flags[id].nil?
            seq = line[7..-1] # "  AAAA----" or "AAA-----" or
            if gap_flags[id].empty? && !id.start_with?('R') # id[0..0] != "R" # seqのほう
              gap_count = seq.rstrip.count("\s")
            end
            seq = seq.tr("\s", "-") # " "を"-"に変更
            gap_flags[id] += seq
          end
        when 2
          break
        end
      end
      return [nil, nil, nil, nil] if sttpos == 0 # 161102に修正

      align_seq = ''
      align_ref = ''
      sttpos   = sttpos - gap_count - 1
      if sttpos < 1
        ref_head = ""
      else
        ref_head = reads.last.last[0..(sttpos-1)]
      end

      seq_head = "-" * sttpos
      gap_flags.each do |flag, seq|
        if flag.start_with?('R') #/^R/ =~ flg
          align_ref = (ref_head + seq)
        else
          align_seq = (seq_head + seq)
        end
      end
      # seqとrefの長さを揃える
      if align_ref.size != align_seq.size
        back_gap = "-" * (align_ref.size - align_seq.size)
        align_seq += back_gap
      end

      # gap_priorityのalignmentのcheck
      check_flg = 0
      if align_ref.count("-") == 0   # deletion check
        check_flg += 1 if /(^\-+)(\w+)(\-+)(\w+)(\-+)$/ =~ align_seq
      elsif /(^\w+)(-+)(\w+)$/ =~ align_ref   # insertion
        check_flg += 1 if /(^\-+)(\w+)(\-+)$/ =~ align_seq
      end

      # ------------------------------------------
      # gap_priorityのalignmentがおかしいときnon_gap_priorityを考慮する
      # ------------------------------------------
      if check_flg == 0
        sttpos  = 0
        gap_count = 0
        count = 0
        non_gap_priority_alignment.split("\n").each do |line|
          count += 1 if line.start_with?('>') # /^>/ =~ line
          case count
          when 0
            next
          when 1
            sttpos = $1.to_i if /^global\/local score.+:(.+)-\d+\)$/ =~ line
            id = line.split("\s").first
            unless non_gap_flags[id].nil?
              seq = line[7..-1] # "  AAAA----" or "AAA-----" or
              if non_gap_flags[id].empty? && !id.start_with?('R') #id[0..0] != "R" # seqのほう
                gap_count = seq.rstrip.count("\s")
              end
              seq = seq.tr("\s", "-") # " "を"-"に変更
              non_gap_flags[id] += seq
            end
          when 2
            break
          end
        end

        align_seq = ''
        align_ref = ''
        sttpos = sttpos - gap_count - 1
        if sttpos < 1
          ref_head = ""
        else
          ref_head = reads.last.last[0..(sttpos-1)]
        end
        seq_head = "-" * sttpos

        non_gap_flags.each do |flag, seq|
          if flag.start_with?('R') #/^R/ =~ flg
            align_ref = (ref_head + seq)
          else
            align_seq = (seq_head + seq)
          end
        end
        # seqとrefの長さを揃える
        if align_ref.size != align_seq.size
          back_gap = "-" * (align_ref.size - align_seq.size)
          align_seq += back_gap
        end
      end

      # ---------------------------------------------------------------------------------------
      #
      # makeing a consensus seq
      #
      # ---------------------------------------------------------------------------------------
      aln = Bio::Alignment.new([align_seq, align_ref])
      consensus = aln.consensus_string(percent_id, gap_mode: -1) # threshold =%id
      consensus, _ = trim_consensus(consensus)
      align_reads_names = []
      align_reads_names << "#{align_seq}\t#{flag_name["seq"]}" << "#{align_ref}\t#{flag_name["ref"]}" << "#{consensus}\tCONSENSUS"

      if trim_flag.nil?
        trim_flag = trim_consensus(consensus)
      else
        trim_consensus_with_flag(consensus, trim_flag)
      end
      return consensus, align_reads_names, count, trim_flag
    end

    private

    def exec_glsearch(reads)
      flags1 = {}
      flags2 = {}
      flag_name = {}

      tmp_seq = Tempfile.new("gls_seq", @temp_path)
      tmp_ref = Tempfile.new("gls_ref", @temp_path)
      reads.each do |read_inf|
        read_flag, sttpos, endpos, seq = read_inf
        flags1["#{read_flag}_#{sttpos}_#{endpos}"[0..5]] = ""
        flags2["#{read_flag}_#{sttpos}_#{endpos}"[0..5]] = ""
        case read_flag
        when :R
          tmp_ref.puts ">#{read_flag}_#{sttpos}_#{endpos}"
          tmp_ref.puts seq.upcase
          flag_name["ref"] = "#{read_flag}_#{sttpos}_#{endpos}"
        else
          tmp_seq.puts ">#{read_flag}_#{sttpos}_#{endpos}"
          tmp_seq.puts seq.upcase
          flag_name["seq"] = "#{read_flag}_#{sttpos}_#{endpos}"
        end
      end
      tmp_seq.flush
      tmp_ref.flush
      cmd = "#{@glsearch_bin} -s #{@glsearch_mat} -g0 -f20 #{tmp_seq.path} #{tmp_ref.path}"
      gap_priority = `${cmd}` # gap優先alignment
      report_error($?, cmd, [tmp_seq, tmp_ref]) unless $?.success?
      cmd = "#{@glsearch_bin} -s #{@glsearch_mat} -f20 #{tmp_seq.path} #{tmp_ref.path}"
      non_gap_priority = `#{cmd}` # gap非優先alignment
      report_error($?, cmd, [tmp_seq, tmp_ref]) unless $?.success?
      tmp_seq.close(true)
      tmp_ref.close(true)
      return { gap_priority_alignment: gap_priority,
               non_gap_priority_alignment: non_gap_priority,
               gap_flags: flags1,
               non_gap_flags: flags2,
               flag_name: flag_name }
    end
  end
end
