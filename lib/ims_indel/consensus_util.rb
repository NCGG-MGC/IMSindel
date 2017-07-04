module IMSIndel::ConsensusUtil
  def trim_first_and_last_3bases(consensus)
    consensus[3..-4] # triming of first 3 base and last 3 base
  end

  def trim_first_3bases(consensus)
    consensus[3..-1] # triming of first 3
  end

  def trim_last_3bases(consensus)
    consensus = consensus[0..-4] # triming of last 3
  end

  def trim_consensus(consensus_str)
    trim_flag = :trim_no
    first3 = consensus_str[0..2].include?("?")
    last3 = consensus_str[-3..-1].include?("?")

    if first3 && last3
      trim_flag = :trim_both
      consensus_str = trim_first_and_last_3bases(consensus_str)
    elsif first3
      trim_flag = :trim_before
      consensus_str = trim_first_3bases(consensus_str)
    elsif last3
      trim_flag = :trim_after
      consensus_str = trim_last_3bases(consensus_str)
    end
    return consensus_str, trim_flag
  end

  def trim_consensus_with_flag(consensus_str, trim_flag)
    case trim_flag
    when :trim_both
      return trim_first_and_last_3bases(consensus_str)
    when :trim_before
      return trim_first_3bases(consensus_str)
    when :trim_after
      return trim_last_3bases(consensus_str)
    end
    return consensus_str
  end
end
