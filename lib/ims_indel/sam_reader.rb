require_relative 'sam_row'

module IMSIndel
  class SAMReader
    PROPER_PAIR    = 0x02
    UNMAPPED       = 0x04
    MATE_UNMAPPED  = 0x08
    REVERSE_STRAND = 0x10
    DUPLICATE      = 0x400

    def initialize(samtools, bam)
      @samtools = samtools
      @bam = bam
    end

    def each(chr:, filter_flag: nil, output_flag: nil, &block)
      command =[@samtools, "view"]
      command << '-F' << filter_flag.to_s if filter_flag
      command << '-f' << output_flag.to_s if output_flag
      command << @bam << chr
      puts command.join(' ')
      IO.popen(command) do |input|
        input.each_line do |line|
          block.call(SAMRow.new(line))
        end
      end
    end
  end
end
