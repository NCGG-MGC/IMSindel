module IMSIndel
  class SAMRow
    attr_accessor :read_name, :flag, :chr, :chrpos, :map_score,
      :map_status, :read_seq, :base_qual, :insert_size

    def initialize(line)
      @read_name, @flag, @chr, chrpos, map_score, @map_status, _, _,
        insert_size, @read_seq, @base_qual,  = line.split("\t", 12)

      @map_score   = map_score.to_i
      @insert_size = insert_size.to_i
      @chrpos   = chrpos.to_i
      @read_seq.upcase!
    end 

    def clip_seq(min_length, min_base_quality)
      if m = /\A(\d+)M(\d+)S\z/.match(@map_status) # 80M20S
        backward_flagment_match(m, min_length, min_base_quality)
      elsif m = /\A(\d+)M(\d+)I(\d+)M(\d+)S\z/.match(@map_status)
        backward_flagment_short_ins(m, min_length, min_base_quality)
      elsif m = /\A(\d+)M(\d+)D(\d+)M(\d+)S\z/.match(@map_status)
        backward_flagment_short_del(m, min_length, min_base_quality)
      elsif m = /\A(\d+)S(?:\d+)M\z/.match(@map_status)
        forward_flagment(m, min_length, min_base_quality)
      elsif m = /\A(\d+)S(?:\d+)M(?:\d+)I(?:\d+)M\z/.match(@map_status)
        forward_flagment(m, min_length, min_base_quality)
      elsif m = /\A(\d+)S(?:\d+)M(?:\d+)D(?:\d+)M\z/.match(@map_status)
        forward_flagment(m, min_length, min_base_quality)

      # just short INDEL
      elsif m = /\A(\d+)M(\d+)I(\d+)M\z/.match(@map_status)
        short_insertion(m, min_length, min_base_quality)
      elsif m = /\A(\d+)M(\d+)D(\d+)M\z/.match(@map_status)
        short_deletion(m, min_length, min_base_quality)
      end
    end

    def acceptable_base_qaul?(base_qual, min_bq)
      average_bq(base_qual) > min_bq
    end

    def average_bq(bq)
      bq.each_byte.inject(0) { |res, base| res += (base - 33) } / bq.size.to_f
    end

    def flagged?(flag)
      @flag.to_i & flag
    end

    private

    def backward_flagment_match(match, min_len, min_bq) # B fragment (match + clip seq)
      len_s = match[2].to_i
      return if len_s <= min_len
      return unless acceptable_base_qaul?(@base_qual[-len_s, len_s], min_bq)

      len_m = match[1].to_i
      [:B, @read_seq[-len_s, len_s], @chrpos + len_m - 1, len_s]
    end

    def backward_flagment_short_ins(match, min_len, min_bq) # B fragment (short INS + clip seq)
      len_s = match[4].to_i
      return if len_s <= min_len
      return unless acceptable_base_qaul?(@base_qual[-len_s, len_s], min_bq)

      len_m1  = match[1].to_i
      len_m2  = match[3].to_i
      [:B, @read_seq[-len_s, len_s], @chrpos + len_m1 + len_m2 - 1, len_s]
    end

    def backward_flagment_short_del(match, min_len, min_bq) # B fragment (short DEL + clip seq)
      len_s = match[4].to_i
      return if len_s <= min_len
      return unless acceptable_base_qaul?(@base_qual[-len_s, len_s], min_bq)

      len_m1  = match[1].to_i
      len_d   = match[2].to_i
      len_m2  = match[3].to_i
      [:B, @read_seq[-len_s, len_s], @chrpos + len_m1 + len_d + len_m2 - 1, nil]
    end

    def forward_flagment(match, min_len, min_bq) # F fragment (clip seq + match)
      len_s = match[1].to_i
      return if len_s <= min_len
      return unless acceptable_base_qaul?(@base_qual[0, len_s], min_bq)

      [:F, @read_seq[0, len_s], @chrpos, len_s]
    end

    def short_insertion(match, min_len, min_bq)
      len_m1  = match[1].to_i
      len_ins = match[2].to_i

      start_pos = len_m1
      end_pos = len_m1 + len_ins - 1
      [:SI, @read_seq[start_pos..end_pos], @chrpos + len_m1 - 1, nil]
    end

    def short_deletion(match, min_len, min_bq)
      len_m1  = match[1].to_i
      len_del = match[2].to_i

      start_pos = len_m1
      end_pos = len_m1 + len_del - 1
      [:SD, @read_seq[start_pos..end_pos], @chrpos + len_m1 - 1, nil]
    end
  end
end
