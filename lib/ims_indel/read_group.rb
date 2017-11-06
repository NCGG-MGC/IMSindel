class IMSIndel::ReadGroup
  attr_reader :backward_clips, :forward_clips, :non_clips

  def initialize(reads, within, support_reads, support_clip_length)
    @backward_clips, backward_non_clips =
       classify(reads.backward_clips + reads.non_clips, within, support_reads, support_clip_length)
    @forward_clips, forward_non_clips =
       classify(reads.forward_clips + reads.non_clips, within, support_reads, support_clip_length)
    @non_clips = make_unique_clips(backward_non_clips, forward_non_clips)
  end

  private

  # grouping reads within 3bp
  def classify(reads, within, support_reads, support_clip_length)
    read_groups = grouping_reads(reads, within, support_reads, support_reads)
    clips = [] #  without SID reads
    non_clips = []

    read_groups.each do |read_group|
      if read_group.all? { |read| read.type == :B || read.type == :F }
        clips << read_group.sort
      else
        non_clips << read_group.sort
      end
    end
    non_clips.uniq!
    return clips, non_clips
  end

  def grouping_reads(reads, within, support_reads, support_clip_length)
    # 1. grouping the reads within 3bp
    sorted_reads = reads.sort{ |a, b| a.start_pos <=> b.start_pos }
    mixed_groups = sorted_reads.inject([]) do |acc, read|
      if acc.empty?
        acc << [read]
      else
        pre_start_pos = acc.last.last.start_pos
        pre_end_pos = acc.last.last.end_pos
        if pre_start_pos <= read.start_pos && read.start_pos <= pre_end_pos + within # within 3bp
          acc.last << read
        else
          acc << [read]
        end
      end
      acc
    end

    # 2. if SD and SI are mixed, major is selected
    read_groups = []
    mixed_groups.each do |group|
      sum_short_insertion = group.count { |i| i.type == :SI }
      sum_short_deletion  = group.count { |i| i.type == :SD }

      # selection of majority of SD and SI
      if sum_short_deletion == 0 && sum_short_insertion == 0
        read_groups << group
      elsif sum_short_deletion < sum_short_insertion
        read_groups << group.find_all { |read| read.type != :SD }
      else
        read_groups << group.find_all { |read| read.type != :SI }
      end
    end

    # 3. check if the number of support reads is more than support_reads
    #    and check any read has soft clip length > support_clip_length
    return read_groups.find_all { |group| group.size >= support_reads }
    #return read_groups.find_all { |group|
      #group.size >= support_reads &&
        #(group.all? { |read| read.soft_clip_len.nil? } || # no sof clip
        #group.any?{ |read| read.soft_clip_len.to_i > support_clip_length })
    #}
  end

  def make_unique_clips(backward_non_clips, forward_non_clips)
    non_clips = []
    used_findex = []

    backward_non_clips.each do |backward_reads|
      tmp = []
      forward_non_clips.each_with_index do |forward_reads, index|
        unless (backward_reads & forward_reads).empty?
          tmp = backward_reads
          forward_reads.each do |ary|
            tmp << ary unless tmp.include?(ary)
          end
          used_findex << index
          break
        end
      end
      if tmp.empty?
        non_clips << backward_reads
      else
        non_clips << tmp
      end
    end
    forward_non_clips.each_with_index do |forward, index|
      non_clips << forward unless used_findex.include?(index)
    end
    return non_clips
  end
end
