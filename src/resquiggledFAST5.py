from Bio import SeqIO
import numpy as np
import os
import io
import h5py
from queue import Queue

import signal_pb2
import nanopolish_pb2
import raw_current_pb2

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
        self._events = None
        self._nanopolish_events = None
        self._raw_current = None
        self._discrete_signal = None
        self._continous_signal = None
        self._add_conversion_data()


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
    

    def get_signal_discrete(self):
        """Returns all discrete signal values associated with the current file.
      
        Parameters
        ---------- 

        Returns
        -------
        signal : np.array
            Discrete signal signal values from file 
        """
        if self._discrete_signal is not None:
            return self._discrete_signal

        signal = self._key_dict_flat['Signal']
    
        self._discrete_signal = np.array(signal)

        return self._discrete_signal

        
    def get_signal_continuos(self):
        """Returns all continuous/raw signal values associated with the current file.
      
        Parameters
        ---------- 

        Returns
        -------
        signal : np.array
            Continuous/ras signal values from file 
        """
        if self._continous_signal is not None:
            return self._continous_signal

        discrete_signal = self.get_signal_discrete()

        info = self.get_general_info()
        rng = info.range
        digitisation = info.digitisation
        raw_unit = rng / digitisation
        offset = info.offset

        self._continous_signal = np.array(list(map(lambda x: (x + offset) * raw_unit, discrete_signal)))
    
        return self._continous_signal
        

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
        events = self.get_events()
        indices = []
        for event in events:
            if event.base == nucleotide:
                indices += [i for i in range(event.start, event.start + event.length)]

        return np.array(indices)


    #TODO
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
            Numpy array containing the indices.
        """
        events = self.get_events()
        indices = []
        for event in events:
            if event.base == nucleotide:
                indices.append((event.start, event.start + event.length))

        return np.array(indices)



    #TODO
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
    

    def get_nucleotide_mean_stdev(self, nuclotide):
        """Returns mean and standard deviation of the specified nucleotide in the signal.
        
        Parameters
        ----------
        nucleotide : string
            Nucleotide representation


        Returns
        -------
        params : tuple
            Tuple of format (mean, stdev)
        """
        events = self.get_events()
        
        mean_sum = 0.0
        stdev_sum = 0.0
        counter = 0

        for event in events:
            if event.base == nucleotide:
                mean_sum  += event.norm_mean
                stdev_sum += event.norm_stdev
                counter   += 1

        return (mean_sum / counter, stdev_sum / counter)


    #TODO
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
        events = self.get_events()
        
        duration_sum = 0
        counter = 0

        for event in events:
            if event.base == nucleotide:
                duration_sum += event.length
                counter      += 1

        return float(duration_sum) / counter


    #TODO
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

    @staticmethod
    def _create_raw_current(signal, fast5_info):
        """MOVE TO SEPARATE FILE"""
        curr = raw_current_pb2.Current()
        
        curr.digitisation = fast5_info.digitisation
        curr.offset = fast5_info.offset
        curr.range = fast5_info.range
        curr.sampling_rate = fast5_info.sampling_rate
        curr.currents.extend(signal)

        return curr


    @staticmethod
    def _create_Event(norm_mean, norm_stdev, start, length, base):
        """MOVE TO SEPARATE FILE"""
        event = signal_pb2.Event()
       
        event.norm_mean  = norm_mean
        event.norm_stdev = norm_stdev
        event.start      = start
        event.length     = length
        event.base       = base
        
        return event
    
    @staticmethod
    def _create_nanopolish_event(events, index, resquiggle_info, fast5_info, aligns):
        """MOVE TO SEPARATE FILE"""
        event_nano = nanopolish_pb2.EventAlign.Event()
        event = events[5]
        increasing = (resquiggle_info.mapped_start < resquiggle_info.mapped_end)

        #could also be clipped_bases_start, end
        event_nano.index = (resquiggle_info.mapped_start + index) if increasing  else (resquiggle_info.mapped_start - index)
        event_nano.level_mean = np.average(list(map(lambda event: np.average(list(map(lambda x: x, event.samples))), events)))
        event_nano.stdv = np.average(list(map(lambda event: np.std(list(map(lambda x: x, event.samples))), events)))
        event_nano.length = np.average(list(map(lambda x: x.length / fast5_info.sampling_rate, events)))
        # event_nano.level_mean = np.average(list(map(lambda x: x, event.samples)))
        # event_nano.stdv = np.std(list(map(lambda x: x, event.samples)))
        # event_nano.length = event.length / fast5_info.sampling_rate
        event_nano.start_idx = resquiggle_info.mapped_start + event.start
        event_nano.end_idx = resquiggle_info.mapped_start + event.start + event.length
        event_nano.samples.extend(event.samples)
        event_nano.standardized_level = event.norm_mean

        kmer = "".join(list(map(lambda x: x.base, events[0:6])))
        if kmer not in aligns:
            event_align = nanopolish_pb2.EventAlign()
            event_align.contig = resquiggle_info.mapped_chrom
            event_align.position = index
            event_align.read_index = 0 #TODO
            event_align.reference_kmer = kmer
            event_align.strand = resquiggle_info.mapped_strand == '+'
            aligns[kmer] = event_align
        event_align = aligns[kmer]

        def reverseComp(char):
            if char == 'T':
                return 'A'
            if char == 'A':
                return 'T'
            if char == 'C':
                return 'G'
            if char == 'G':
                return 'C'
        #either model or reference is from Analyses.BaseCall
        if 'c' in resquiggle_info.mapped_chrom:
            event_align.model_kmer = "".join(list(map(reverseComp, kmer[::-1])))
        else:
            event_align.model_kmer = kmer
        event_align.events.extend([event_nano])
        event_align.model_mean = np.average(list(map(lambda x: x.level_mean, event_align.events)))
        event_align.model_stdv = np.average(list(map(lambda x: x.stdv, event_align.events)))

        return event_align
    @staticmethod
    def _create_ResquiggleInfo(attrs):
        """MOVE TO SEPARATE FILE"""
        info = signal_pb2.ResquiggleInfo()
        
        info.clipped_bases_end   = int(attrs['clipped_bases_end'])
        info.clipped_bases_start = int(attrs['clipped_bases_start'])
        info.mapped_chrom        = attrs['mapped_chrom']
        info.mapped_end          = int(attrs['mapped_end'])
        info.mapped_start        = int(attrs['mapped_start'])
        info.mapped_strand       = attrs['mapped_strand']
        info.num_deletions       = int(attrs['num_deletions'])
        info.num_insertions      = int(attrs['num_insertions'])
        info.num_matches         = int(attrs['num_matches'])
        info.num_mismatches      = int(attrs['num_mismatches'])
       
        return info

    @staticmethod
    def _create_Fast5Info(template_attrs, channel_id_attrs):
        """MOVE TO SEPARATE FILE"""
        info = signal_pb2.Fast5Info()
       
        info.lower_lim          = float(template_attrs['lower_lim'])
        info.norm_type          = template_attrs['norm_type']
        info.outlier_threshold = float(template_attrs['outlier_threshold'])
        info.rna                = bool(template_attrs['rna'])
        info.scale              = float(template_attrs['scale'])
        info.shift              = float(template_attrs['shift'])
        info.signal_match_score = float(template_attrs['signal_match_score'])
        info.status             = template_attrs['status']
        info.upper_lim          = float(template_attrs['upper_lim'])

        info.channel_number     = channel_id_attrs['channel_number']
        info.digitisation       = float(channel_id_attrs['digitisation'])
        info.offset             = float(channel_id_attrs['offset'])
        info.range              = float(channel_id_attrs['range'])
        info.sampling_rate      = float(channel_id_attrs['sampling_rate'])

        return info

    @staticmethod
    def _serialize_Message(message, path):
        """MOVE TO SEPARATE FILE"""
        f = open(path, "wb")
        f.write(message.serializeToString())
        f.close()
        return 


    def get_resquiggle_info(self):
        """Returns resquiggle information as ResquiggleInfo object.

        Parameters
        ----------

        Return
        info : ResquiggleInfo
            Object containing the information as object variables
        """
        attrs = self._key_dict_flat['Alignment'].attrs
        return ResquiggledFAST5._create_ResquiggleInfo(attrs)


    def get_general_info(self):
        """Returns general information about the sequence as well as tombo parameters 
        as a Fast5Info object.

        Parameters
        ----------
        
        Returns
        -------
        info : Fast5Info

        """
        template_attrs = self._key_dict_flat['BaseCalled_template'].attrs
        channel_id_attrs = self._key_dict_flat['channel_id'].attrs
        return ResquiggledFAST5._create_Fast5Info(template_attrs, channel_id_attrs)

    def get_nanopolish_events(self):
        """Returns all Event objects in sequential order associated with the fast5 file converted to nanopolish output.
        
        Parameters
        ----------

        Returns
        -------
        events : [Event] ili numpy.array
            All events in sequential order

        """
        if self._nanopolish_events is not None:
            return self._nanopolish_events

        events = self.get_events()
        nano_events = []
        aligns = {}
        for i in range(3,len(events)-3):
            nano_events.append(self._create_nanopolish_event(events[i-3:i+3], i-3, self.get_resquiggle_info(), self.get_general_info(), aligns))

        self._nanopolish_events = nano_events
        return self._nanopolish_events

    def get_raw_current_data(self):
        """Returns Current object retrieved from the signal in the Fast5 file
        
        Parameters
        ----------

        Returns
        -------
        events : Current
            Current file 

        """
        if self._raw_current is not None:
            return self._raw_current

        signal = self.get_signal_continuos()

        self._raw_current = ResquiggledFAST5._create_raw_current(signal, self.get_general_info())
        return self._raw_current
        

    def get_events(self):
        """Returns all Event objects in sequential order associated with the fast5 file.
        
        Parameters
        ----------

        Returns
        -------
        events : [Event] ili numpy.array
            All events in sequential order

        """
        if self._events is not None:
            return self._events

        alignment = self._key_dict_flat['BaseCalled_template']
        events = np.array(alignment.get('Events'))
        signal = self.get_signal_continuos()                            # events is a vector of shape (x, ) 

        events_list = []
        for row in events:
            event = ResquiggledFAST5._create_Event(*row)
            event.samples.extend(signal[event.start:(event.start+event.length)])
            events_list.append(event)
            
    
        self._events = np.array(events_list)

        return self._events


    def _convert_to_raw(self, signals):
        """Converts given np.array of discrete (int) signals in to their continuous/raw (float) counterparts.

        Parameters
        ----------
        signals : np.array
            Array of discrete signals
        
        Returns
        -------
        raw_signals : np.array
            Array of raw/continuous signals
        """
        return np.vectorize(self._convert_single_to_raw, signals)


    def _convert_single_to_raw(self, disc_signal):
        """Helper function which converts a discrete signal into a continuous/raw one.
        The transformation is given by the following formula:
            (disc_signal + offset) / raw_unit
        where raw_unit is computed as:
            raw_unit = range / digitisation
        
        Offset, range and digitisation are all parameters specified in the fast5 file.
        
        Parameters
        ----------
        disc_signal : int
            Discrete representation of signal value

        Returns
        -------
        cont_signal : float
            Continuous representation of signal value
        """
        return (disc_signal + self._offset) * self._raw_unit 


    def _add_conversion_data(self):
        """Reads variables needed for signal conversion (digitisation, offset and range) into instance variables.
        
        Parameters
        ----------

        Returns
        -------
        """
        attrs = self._key_dict_flat['channel_id'].attrs
        self._digitisation = float(attrs['digitisation'])
        self._offset       = float(attrs['offset'])
        self._range        = float(attrs['range'])
        self._raw_unit     = self._range / self._digitisation
        #self._raw_unit = self._digitisation / self._range

    def _generate_key_dict(self):
        """Generates key dictionary for accessing all groups in fast5 file for easier access in other methods.
        Two dictionaries are generated:
            _key_dict_flat : (key string) -> (HDF5 group object)
            _key_dict_hierarchical : (HDF5 group object) -> ([key string])

        Parameters
        ----------

        Returns 
        -------
        """
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
                
        #"""UNCOMMENT TO VIEW DICT ELEMENTS
        print("Hierarchical dict representation")        
        for key, value in self._key_dict_hierarchical.items():
            print(key, value)
        print(" ")
        print("Flat dict representation: ")
        
        for key, value in self._key_dict_flat.items():
            print(key, value)
        #"""

        return


    def _extract_file_format(self, file_format):
        """Extracts the specified file format from fast5 file. If no file with the specified format is found,
        an Error is raised.

        Parameters
        ----------
        file_format : string
            File format as string.

        Returns
        -------
        seq : Bio.Seq.Seq
            Sequence object representation of the specified file format.


        Raises
        ------
        KeyError
            Raised if there is no file with the specified format
        """
        result = None
        
        try:
            dataset = self._key_dict_flat[file_format]
            string = dataset[()]
            result = list(SeqIO.parse(io.StringIO(string), file_format.lower()))    #for simple conversion
                                                                                    #we can remove list 

        except KeyError as key_error:
            print("This fast5 file does not contain a {} file.".format(file_format.lower()))
        
        return result
