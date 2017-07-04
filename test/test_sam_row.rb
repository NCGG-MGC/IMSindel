require 'test/unit'
require_relative '../lib/ims_indel/sam_row'

class TestBackwardMatch < Test::Unit::TestCase
  def setup
    @sam = IMSIndel::SAMRow.new("QNAME\t2\t22\t10\t50\t3M4S\t=\t0\t10\tATGCGTA\tABCDEFG")
  end

  def test_clip_seq
    type, read, pos = @sam.clip_seq(1, 1)
    assert_equal :B, type
    assert_equal 'CGTA', read
    assert_equal 12, pos
  end

  def test_clip_seq_min_length_good
    assert_not_equal nil, @sam.clip_seq(3, 1)
  end
  def test_clip_seq_min_length_bad
    assert_equal nil, @sam.clip_seq(4, 1)
  end
  def test_clip_seq_base_qual_bood
    assert_not_equal nil, @sam.clip_seq(3, 36.4)
  end
  def test_clip_seq_base_qual_bad
    assert_equal nil, @sam.clip_seq(3, 36.5)
  end
end

class TestBackwardShortIns < Test::Unit::TestCase
  def setup
    @sam = IMSIndel::SAMRow.new("QNAME\t2\t22\t10\t50\t1M2I3M4S\t=\t0\t10\tATGCGTA\tABCDEFG")
  end

  def test_clip_seq
    type, read, pos = @sam.clip_seq(1, 1)
    assert_equal :B, type
    assert_equal 'CGTA', read
    assert_equal 13, pos
  end

  def test_clip_seq_min_length_good
    assert_not_equal nil, @sam.clip_seq(3, 1)
  end
  def test_clip_seq_min_length_bad
    assert_equal nil, @sam.clip_seq(4, 1)
  end
  def test_clip_seq_base_qual_bood
    assert_not_equal nil, @sam.clip_seq(3, 36.4)
  end
  def test_clip_seq_base_qual_bad
    assert_equal nil, @sam.clip_seq(3, 36.5)
  end
end

class TestBackwardShortDel < Test::Unit::TestCase
  def setup
    @sam = IMSIndel::SAMRow.new("QNAME\t2\t22\t10\t50\t1M2D3M4S\t=\t0\t10\tATGCGTA\tABCDEFG")
  end

  def test_clip_seq
    type, read, pos = @sam.clip_seq(1, 1)
    assert_equal :B, type
    assert_equal 'CGTA', read
    assert_equal 15, pos
  end

  def test_clip_seq_min_length_good
    assert_not_equal nil, @sam.clip_seq(3, 1)
  end
  def test_clip_seq_min_length_bad
    assert_equal nil, @sam.clip_seq(4, 1)
  end
  def test_clip_seq_base_qual_bood
    assert_not_equal nil, @sam.clip_seq(3, 36.4)
  end
  def test_clip_seq_base_qual_bad
    assert_equal nil, @sam.clip_seq(3, 36.5)
  end
end

class TestForwardMatch < Test::Unit::TestCase
  def setup
    @sam = IMSIndel::SAMRow.new("QNAME\t2\t22\t10\t50\t4S3M\t=\t0\t10\tATGCGTA\tDEFGABC")
  end

  def test_clip_seq
    type, read, pos = @sam.clip_seq(1, 1)
    assert_equal :F, type
    assert_equal 'ATGC', read
    assert_equal 10, pos
  end

  def test_clip_seq_min_length_good
    assert_not_equal nil, @sam.clip_seq(3, 1)
  end
  def test_clip_seq_min_length_bad
    assert_equal nil, @sam.clip_seq(4, 1)
  end
  def test_clip_seq_base_qual_bood
    assert_not_equal nil, @sam.clip_seq(3, 36.4)
  end
  def test_clip_seq_base_qual_bad
    assert_equal nil, @sam.clip_seq(3, 36.5)
  end
end

class TestForwardShortIns < Test::Unit::TestCase
  def setup
    @sam = IMSIndel::SAMRow.new("QNAME\t2\t22\t10\t50\t4S3M2I1M\t=\t0\t10\tATGCGTA\tDEFGABC")
  end

  def test_clip_seq
    type, read, pos = @sam.clip_seq(1, 1)
    assert_equal :F, type
    assert_equal 'ATGC', read
    assert_equal 10, pos
  end

  def test_clip_seq_min_length_good
    assert_not_equal nil, @sam.clip_seq(3, 1)
  end
  def test_clip_seq_min_length_bad
    assert_equal nil, @sam.clip_seq(4, 1)
  end
  def test_clip_seq_base_qual_bood
    assert_not_equal nil, @sam.clip_seq(3, 36.4)
  end
  def test_clip_seq_base_qual_bad
    assert_equal nil, @sam.clip_seq(3, 36.5)
  end
end

class TestForwardShortDel < Test::Unit::TestCase
  def setup
    @sam = IMSIndel::SAMRow.new("QNAME\t2\t22\t10\t50\t4S3M2D1M\t=\t0\t10\tATGCGTA\tDEFGABC")
  end

  def test_clip_seq
    type, read, pos = @sam.clip_seq(1, 1)
    assert_equal :F, type
    assert_equal 'ATGC', read
    assert_equal 10, pos
  end

  def test_clip_seq_min_length_good
    assert_not_equal nil, @sam.clip_seq(3, 1)
  end
  def test_clip_seq_min_length_bad
    assert_equal nil, @sam.clip_seq(4, 1)
  end
  def test_clip_seq_base_qual_bood
    assert_not_equal nil, @sam.clip_seq(3, 36.4)
  end
  def test_clip_seq_base_qual_bad
    assert_equal nil, @sam.clip_seq(3, 36.5)
  end
end

class TestShortInsertion < Test::Unit::TestCase
  def setup
    @sam = IMSIndel::SAMRow.new("QNAME\t2\t22\t10\t50\t3M2I1M\t=\t0\t10\tATGCGTA\tDEFGABC")
  end

  def test_clip_seq
    type, read, pos = @sam.clip_seq(1, 1)
    assert_equal :SI, type
    assert_equal 'CG', read
    assert_equal 12, pos
  end
end

class TestShortDeletion < Test::Unit::TestCase
  def setup
    @sam = IMSIndel::SAMRow.new("QNAME\t2\t22\t10\t50\t3M2D1M\t=\t0\t10\tATGCGTA\tDEFGABC")
  end

  def test_clip_seq
    type, read, pos = @sam.clip_seq(1, 1)
    assert_equal :SD, type
    assert_equal 'CG', read
    assert_equal 12, pos
  end
end
