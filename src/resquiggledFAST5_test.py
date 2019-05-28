import resquiggledFAST5 as rs
import sys
if __name__ == "__main__":
    filename = "./data/dac.fast5"
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    f = rs.ResquiggledFAST5(filename)
    print(f.get_nucleotide_positions.__doc__)
    #print(f.get_fastq())
    print(f.get_events())
    print(f.get_signal_discrete())
    print(f.get_signal_continuos())
    print(f.get_general_info())
    print(f.get_nucleotide_intervals('T'))
    print(f.get_nanopolish_events())
