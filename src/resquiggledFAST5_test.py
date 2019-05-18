import resquiggledFAST5 as rs

if __name__ == "__main__":
    f = rs.ResquiggledFAST5("./data/db6b45aa-5d21-45cf-a435-05fb8f12e839.fast5")
    print(f.get_nucleotide_positions.__doc__)
    #print(f.get_fastq())
    print(f.get_events())
    print(f.get_signal_discrete())
    print(f.get_signal_continuos())
    print(f.get_general_info())
    print(f.get_nucleotide_intervals('T'))
