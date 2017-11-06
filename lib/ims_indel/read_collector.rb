require_relative 'sam_reader'
require_relative 'read'

module IMSIndel
  class ReadCollector

    attr_reader :backward_clips,        # back soft clipありのlong indel用(readがread_seq)
                :forward_clips,         # forward soft clipありのlong indel用(readがread_seq)
                :non_clips,             # soft clipなしのshort indel用
                :unmapped_reads,
                :pos2max_backward_clip, # back soft clipありのlong indel用(readがclip_seq) { clip_pos => clipseq(最長のもの),
                :pos2max_forward_clip,  # forward soft clipありのlong indel用(readがclip_seq) { clip_pos => clipseq(最長のもの),
                :avg_insert_size,
                :sd,
                :max_read_size,         # 最大完全matched read長
                :min_insert_size,
                :max_insert_size

    def initialize(samtools)
      @samtools = samtools
    end

    def collect_clip_indel_reads(bam, chr, map_quality_threshold,
                                 base_quality_threshold, clip_length)

      max_read_size = 0
      non_clips = []
      backward_clips = []
      forward_clips = []
      insert_sizes = []
      pos2max_backward_clip = {}
      pos2max_forward_clip = {}

      sam = SAMReader.new(@samtools, bam)
      sam.each(chr: chr, filter_flag: SAMReader::DUPLICATE,
                         output_flag: SAMReader::PROPER_PAIR) do |line|

        if line.map_score <= map_quality_threshold || line.read_seq.include?('N')
          next
        end

        if line.map_status =~ /\A(\d+)M\z/
           n = $1.to_i
           max_read_size = n if n > max_read_size
        end

        type, clipped_seq, clipped_pos, soft_clip_len = line.clip_seq(clip_length, base_quality_threshold)
        case type
        when nil
          next
        when :SI
          non_clips << Read.new(type: type, start_pos: clipped_pos, end_pos: clipped_pos, seq: line.read_seq, soft_clip_len: soft_clip_len)
        when :SD
          start_pos = clipped_pos
          end_pos = start_pos + clipped_seq.size - 1
          non_clips << Read.new(type: type, start_pos: start_pos, end_pos: end_pos, seq: line.read_seq, soft_clip_len: soft_clip_len)
        when :B
          backward_clips << Read.new(type: type, start_pos: clipped_pos, end_pos: clipped_pos, seq: line.read_seq, soft_clip_len: soft_clip_len)
          if !pos2max_backward_clip.has_key?(clipped_pos) || pos2max_backward_clip[clipped_pos].size < clipped_seq.size
            pos2max_backward_clip[clipped_pos] = clipped_seq
          end
        when :F
          forward_clips << Read.new(type: type, start_pos: clipped_pos, end_pos: clipped_pos, seq: line.read_seq, soft_clip_len: soft_clip_len)
          if !pos2max_forward_clip.has_key?(clipped_pos) || pos2max_forward_clip[clipped_pos].size < clipped_seq.size
            pos2max_forward_clip[clipped_pos] = clipped_seq
          end
        end
        insert_sizes << line.insert_size if line.insert_size > 0
      end

      calc_stats(insert_sizes)
      @max_read_size = max_read_size
      @non_clips = non_clips
      @backward_clips = backward_clips
      @forward_clips = forward_clips
      @pos2max_backward_clip = pos2max_backward_clip
      @pos2max_forward_clip = pos2max_forward_clip

      insert_sizes.empty?
    end

    def collect_unmapped_reads(bam, chr, map_quality_threshold, base_quality_threshold)
      mate_unmapped_read_names = collect_mate_unmapped_read_name(bam, chr, map_quality_threshold)
      puts "mate_unmapped_read_names: #{mate_unmapped_read_names.size}"
      unmapped_reads = []

      sam = SAMReader.new(@samtools, bam)
      sam.each(chr: chr, filter_flag: SAMReader::DUPLICATE,
                         output_flag: SAMReader::UNMAPPED) do |line|

        next unless line.acceptable_base_qaul?(line.base_qual, base_quality_threshold)

        start_pos = 0
        end_pos = 0
        read_seq = line.read_seq
        type = mate_unmapped_read_names[line.read_name]
        if type == :forward
          start_pos = line.chrpos + @max_read_size - @max_insert_size
          end_pos = chrpos
        elsif type == :backward
          start_pos = line.chrpos + @min_insert_size
          end_pos = line.chrpos + @max_insert_size
          read_seq = rev_comp(read_seq)
          next unless read_seq
        end
        #puts [start_pos, end_pos, read_seq]
        unmapped_reads << Read.new(type: :U, start_pos: start_pos, end_pos: end_pos, seq: read_seq)
      end
      unmapped_reads.sort! { |a, b| a.start_pos <=> b.start_pos }
      @unmapped_reads = unmapped_reads
    end

    # reverse complement of sequences
    def rev_comp(seq) # AAAGG-G => CCCTT-T
      # remove N start and last bases
      seq = seq.sub(/^N+/, '')
      seq.sub!(/N+$/, '')
      return nil if seq.include?("N".freeze)
      seq.upcase!
      seq.reverse!
      seq.tr!('ATGC'.freeze, 'TACG'.freeze)
      return seq
    end

    private

    def calc_stats(insert_sizes)
      @avg_insert_size = calc_avg_insert(insert_sizes)
      var_insert_size = calc_var_insert(insert_sizes, @avg_insert_size)
      @sd  = Math.sqrt(var_insert_size)

      # insert size: +/-3SD
      @min_insert_size = (avg_insert_size - 3 * @sd).to_i
      @max_insert_size = (var_insert_size + 3 * @sd).to_i
    end

    # ---------------------------------------------------------------------------------------
    # Avg, Variance
    # ---------------------------------------------------------------------------------------
    def calc_avg_insert(array)
      array.inject(0) { |res, i| res += i } / array.size.to_f
    end

    def calc_var_insert(array, avg)
      array.inject(0) { |res, i| res += (i - avg) ** 2 } / array.size.to_f
    end

    def collect_mate_unmapped_read_name(bam, chr, map_quality_threshold)
      unmapped_read_names = {}
      sam = SAMReader.new(@samtools, bam)
      sam.each(chr: chr, filter_flag: SAMReader::DUPLICATE,
                         output_flag: SAMReader::MATE_UNMAPPED) do |line|
        next if line.map_score <= map_quality_threshold
        unmapped_read_names[line.read_name] = 
          line.flagged?(SAMReader::REVERSE_STRAND) ? :backward : :forward
      end
      unmapped_read_names
    end
  end
end
