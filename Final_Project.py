class final_project():

    #The function transforms the input file to a dictionary where the key is the header and the value is the sequence
    def __init__(self,dna2fasta):
        self.dna2fasta = dna2fasta
        self.dict = {}
        frame_reader = open (self.dna2fasta)
        for line in frame_reader:
            line = line.strip("\n")
            if ">" in line:
                header = line
                self.dict[header] = ""
            else:
                self.dict[header] += line
        frame_reader.close()
    
    #The function is counting number of records in the file
    def records_count(self):
        number_of_records = len(self.dict)
        print("How many records are in the multi-FASTA file: %d \n"% number_of_records) #Question1: How many records are in the multi-FASTA file?
        return number_of_records

    #The function is checking lenght of the each record and adequate header of the record
    def records_lenght(self):
        lenght_dict = {}
        for key,value in self.dict.items():
            lenght_dict[key] = len(value)
        
        seq_lenght = lenght_dict.values()

        max_seq_lenght = max(seq_lenght)
        min_seq_lenght = min(seq_lenght)

        find_max_seq_lenght = [l for l in lenght_dict if lenght_dict[l] == max_seq_lenght]
        find_min_seq_lenght = [l for l in lenght_dict if lenght_dict[l] == min_seq_lenght]

        #Question2 and Question 3: What is the length of the longest/shortest sequence in the file?
        print ("The lenght of the longest senquence is: %d \n"% max_seq_lenght, 
               "The indentifier of the longest sequence is: %d \n"% len(find_max_seq_lenght))
        print ("The lenght of the shortest senquence is: %d \n"% min_seq_lenght,
                "The indentifier of the shortest sequence is: %d \n"%len(find_min_seq_lenght))
        return lenght_dict

    #The function is finding start position and length for each read frame for a given sequence (dna: sequence, string)   
    def find_position(self, dna):
       
        start_code = "ATG"
        stop_codes = ["TAA", "TAG", "TGA"]

        pos_dict = {} 
        
        for i in range(3): 
            pos = []
            #generate frames
            if i == 0:
                frame = [dna[j:j+3] for j in range(i, len(dna), 3)]
            else:
                frame = [dna[:i]] + [dna[j:j+3] for j in range(i, len(dna), 3)]
            
            #finding all possible start and stop position
            start_pos = []
            stop_pos = []
            try:
                index_start_pos = [m for m, y in enumerate(frame) if y == start_code] 
                start_pos += index_start_pos #start position is "ATG"
            except ValueError:
                pos.append((-1, 0)) #if no "ATG", return -1 as start
                continue
 
            for stop_code in stop_codes:
                try:
                    
                    index_stop_code = [n for n, x in enumerate(frame) if x == stop_code and n > min(start_pos)]
                    stop_pos += index_stop_code
                except ValueError:
                    continue

            if len(stop_pos) == 0: # add -1 as start position when no stop code find
                 pos.append((-1, 0))
            else:
                #find the closest paired code
                 while len(start_pos) != 0:
                     start = min(start_pos)
                     try:
                         end = min([stop for stop in stop_pos if stop > start])
                     except ValueError:
                         break     
                    #adding start position and lenght
                     s_pos = len("".join(frame[:start])) + 1
                     pos.append((s_pos, (end - start + 1)*3))
                     start_pos.remove(start) 
            pos_dict["frame%d"%(i+1)] = pos 
        return pos_dict
    
    # The function return all the ORF information with start posiotion and length of ORF in read frame 1, 2 and 3
    # the values for frame 1, 2, and are represented as pairs of tuple in a list
    # for example: {"header1": {"frame1":[(0, 100)], "frame2":[(20, 400)], "frame3":[(-1, 0)]},...} 
    # represents for header1:
    # frame1- start position: 0, length of orf: 100
    # frame2- start position: 20, length of orf: 400
    # frame3- No ORF detected
    def orf_identifier (self):
        orf = {}
        for header, dna_seq in self.dict.items(): 
            pos = self.find_position(dna_seq)
            orf[header] = pos
        # finding the header for Question7:
        id_key = [key for key in orf if "gi|142022655|gb|EQ086233.1|16" in key]
        idx = id_key[0]

        # finding the list of the frames for Question4-7:  
        frame1, frame2, frame3, all_frames, id_frames = [], [], [], [], []
        for key, dict_value in orf.items():
            frame1 += dict_value["frame1"]
            frame2 += dict_value["frame2"]
            frame3 += dict_value["frame3"]
            frames = dict_value["frame1"] + dict_value["frame2"] + dict_value["frame3"]
            all_frames += frames
            if key == idx:
                id_frames = dict_value["frame1"] + dict_value["frame2"] + dict_value["frame3"]

        #Question4: What is the lenght of the longest ORF appearing in reading frame2?
        frame2_max_length = max(frame2, key = lambda x: x[1])
        print("Q4: The length of longest ORF in frame2: %d\n"%frame2_max_length[1])

        #Question5: What is the starting position of the longest ORF in reading frame3?
        frame1_max_length_pos = max(frame3, key = lambda x: x[1])
        print("Q5: The start position of longest ORF in frame1: %d\n"%frame1_max_length_pos[0])

        #Question6: What is the lenght of the longest ORF appearing in all sequence and reading frames?
        max_length = max(all_frames, key = lambda x: x[1])
        print("Q6: The longest ORF of all frames and sequences: %d\n"%max_length[1])

        #Question7: What is the length of the longest forward ORF that appears in the 
        #sequence with the identifier gi|142022655|gb|EQ086233.1|16?
        max_length_id = max(id_frames, key = lambda x: x[1])
        print("Q7: The length of longest ORF for ", idx, "is: %d \n" %max_length_id[1])

        return orf

    #The function is finding and counting repeats for each dna sequence; dna: sequence, string; n: number of repeats, int
    def repeats_finder(self, dna, n):
        
        repeats = {}
        for i in range(0, len(dna)):
            repeat = dna[i:i+n] 
            if len(repeat) == n:
                if repeat not in repeats:
                    repeats [repeat] = 1 
                else:
                    repeats[repeat] = repeats.get(repeat) + 1
        return repeats

    #The function generates repeats with counts for each record (repeats_set) and for the whole file (combined_repeats); n: number of repeats, int
    def repeats_identifier(self, n):
         
        repeats_set = {}
        for header, dna_seq in self.dict.items():
            repeats = self.repeats_finder(dna_seq, n)
            repeats_set[header] = repeats 
       
        combined_repeats = {}
        for dict_value in repeats_set.values():
            for key in dict_value:
                if key not in combined_repeats:
                    combined_repeats[key] = dict_value[key]
                else:
                    combined_repeats[key] = combined_repeats.get(key) + dict_value[key]
      
        #Question8: Find the most frequently occuring repeat of lenght 6 in all squence and how many times does it occur?            
        if n == 6:
            most_freq_6 = max (combined_repeats.values())
            print("Q8: The most frequently repeats occur: %d times \n"%most_freq_6)
        
        #Question9: How many different 12-base sequences occur Max time?
        if n == 12:
            count_most_freq_12 = len([value for value in combined_repeats.values()\
                             if value == max(combined_repeats.values())])
        
            print("Q9: The number of different 12-base sequences occur max times: %d \n"%count_most_freq_12)

        #Question10: Which one of the following repeats of the lenght 7 has a maximum number of occurrences? 
        if n == 7:
            most_freq_7 = max (combined_repeats.values())
            most_freq_7_seq = [key for key in combined_repeats if \
                       combined_repeats[key] == max(combined_repeats.values())]
            print("Q10: The following repeats occured most frequently: \n"% most_freq_7_seq)

x = final_project("C:/Users/Milica Jokovic/Documents/Course_Python_for_Genomic_Data_Science/dna2.fasta")
lenght = x.records_lenght()
records = x.records_count()
orf = x.orf_identifier()
rep6 = x.repeats_identifier(6)
rep12 = x.repeats_identifier(12)
rep7 = x.repeats_identifier(7)

