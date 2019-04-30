from Bio import SeqIO
import numpy as np
import os
import io
import h5py
from queue import Queue

class ResquiggledFAST5():
    """ <2do> dokumentacija"""

    def __init__(self, path_to_file):
        if not os.path.isfile(path_to_file):
            raise FileNotFoundError("The file {} could not be found".format(path_to_file))

        self._path_to_file = path_to_file
        self._file_handle = h5py.File(path_to_file, 'r')
        self._key_dict_hierarchical = dict()
        self._key_dict_flat = dict()
        self._generate_key_dict()
        self._sequence = None

    def get_fasta(self):
        """Returns fasta reads read within the fast5 file.

        Parameters
        ----------

        Returns
        -------
        fasta : SeqRecord list
            List of SeqRecords objects
        """
        return self._extract_file_format('Fasta')
    
    def get_fastq(self):
        """Returns fastq reads read within the fast5 file.

        Parameters
        ----------

        Returns
        -------
        fastq : SeqRecord list
            List of SeqRecords objects
        """
        return self._extract_file_format('Fastq')
    
    def get_fast5(self):
        """Returns fast5 read values associated with the current file.
      
        Parameters
        ---------- 

        Returns
        -------
        fast5
            Fast5 readings from file ((FORMAT TBD))
        """
       
        pass
    
    def get_nucleotide_positions(self, nucleotide):
        """Returns all indices of positions where the specified nucleotide has been detected.
        
        Parameters
        ----------
        nucleotide : string
            Nucleotide name

        Returns
        -------
        indices : numpy.array
            Numpy vector containing the indices.
        """

        pass

    def get_kmer_positions(self, kmer):
        """Returns all indices of positions where the specified kmer has been detected.
        
        Parameters
        ----------
        kmer : string
            Kmer representation

        Returns
        -------
        indices : numpy.array
            Numpy vector containing the indices.
        """
        pass

    def get_nucleotide_intervals(self, nucleotide):
        """Returns all start and end indices of intervals where the specified nucleotide has been detected.
        
        Parameters
        ----------
        nucleotide : string
            Nucleotide representation


        Returns
        -------
        indices : numpy.array
            Numpy vector containing the indices.
        """
        pass

    def get_kmer_intervals(self, kmer):
        """Returns all start and end indices of intervals where the specified kmer has been detected.
        
        Parameters
        ----------
        kmer : string
            Kmer representation

        Returns
        -------
        indices : numpy.array
            Numpy vector containing the indices.
        """
        pass
    
    def get_nucleotide_mean_sd(self, nuclotide):
        """Returns mean and standard deviation of the specified nucleotide in the signal.
        
        Parameters
        ----------
        nucleotide : string
            Nucleotide representation


        Returns
        -------
        params : tuple
            Tuple of format (mean, std_dev)
        """
        pass

    def get_kmer_mean_sd(self, kmer):
        """Returns mean and standard deviation of the specified kmer in the signal.
        
        Parameters
        ----------
        kmer : string
            Kmer representation


        Returns
        -------
        params : tuple
            Tuple of format (mean, std_dev)
        """
        pass

    def get_nucleotide_average_duration(self, nucleotide):
        """Returns average duration of the specified nucleotide in the signal.
        
        Parameters
        ----------
        nucleotide : string
            Nucleotide representation


        Returns
        -------
        duration : float
            Duration of nucleotide
        """
        pass
    
    def get_kmer_average_duration(self, kmer):
        """Returns average duration of the specified kmer in the signal.
        
        Parameters
        ----------
        kmer : string
            Kmer representation


        Returns
        -------
        duration : float
            Duration of kmer
        """
        pass
    
    def _generate_key_dict(self):
        queue = Queue()
        queue.put(self._file_handle)

        while not queue.empty():
            tmp = queue.get()
            tmp_keys = []
            try:
                for key in tmp.keys():
                    tmp_keys.append(key)
                    queue.put(tmp[key])
                    self._key_dict_flat[key] = tmp[key]

                self._key_dict_hierarchical[tmp] = tmp_keys
            
            except AttributeError as attr_error:
                continue
                
        """UNCOMMENT TO VIEW DICT ELEMENTS
        print("Hierarchical dict representation")        
        for key, value in self._key_dict_hierarchical.items():
            print(key, value)

        print("Flat dict representation: ")
        
        for key, value in self._key_dict_flat.items():
            print(key, value)
        """

        return


    def _extract_file_format(self, file_format):
        result = None
        
        try:
            fastq_dataset = self._key_dict_flat[file_format]
            fastq_string = fastq_dataset[()]
            result = list(SeqIO.parse(io.StringIO(fastq_string), file_format.lower()))

        except KeyError as key_error:
            print("This fast5 file does not contain a {} file.".format(file_format.lower()))
        
        return result