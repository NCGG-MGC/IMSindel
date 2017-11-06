class IMSIndel::Read
  attr_reader :type, :start_pos, :end_pos, :seq, :soft_clip_len, :depth

  def initialize(type:, start_pos:, end_pos:, seq:, soft_clip_len: nil, depth: 0)
    @type = type
    @start_pos = start_pos
    @end_pos = end_pos
    @seq = seq
    @soft_clip_len = soft_clip_len
    @depth = depth
  end

  # 同値判定をちゃんとしておかないと[Read(), ...].uniqあたりでおかしくなる
  def <=>(other)
    c = @type <=> other.type
    return c unless c == 0
    c = @start_pos <=> other.start_pos
    return c unless c == 0
    c = @end_pos <=> other.end_pos
    return c unless c == 0
    c = @seq <=> other.seq
    return c unless c == 0
    return @depth <=> other.depth
  end

  def ==(other)
    eql?(other)
  end

  # hashに入れるためeql? hashをオーバーライド
  def eql?(other)
    return true if equal?(other)
    @type == other.type &&
      @start_pos == other.start_pos &&
      @end_pos == other.end_pos &&
      @seq == other.seq &&
      @depth == other.depth
  end

  def hash
    # もし値を変更するようならキャッシュUPDATE
    @hash_val ||= [@type, @start_pos, @end_pos, @seq, @depth].hash
  end
end
