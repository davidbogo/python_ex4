# Author: David Bogoslavsky
import sys
import math
import random
import re
import os
import csv
import json
import Bio.SeqIO as seqIO
import Bio.Seq as Seq
import matplotlib.pyplot as plt
import numpy as np
import time

#this function check if the dna sequence is actually dna sequence
def check_dna (dna_seq):
    # check if every base is a real dna base
    for i in dna_seq:
        set_dna = {"a", "A", "c", "C", "g", "G", "t", "T"}
        # if not return false
        if i not in set_dna:
            return False
    return True

#this function check if the rna sequence is actually rna sequence
def check_rna (rna_seq):
    # check if every base is a real rna base
    for i in rna_seq:
        set_rna = {"a", "A", "c", "C", "g", "G", "u", "U"}
        #if not return false
        if i not in set_rna:
            return False
    return True

# this class represent a polymerase, dna or rna one.
class Polymerase:
    #initial a polymerase
    def __init__(self, type, error_rate = 0):
        #check if the type is DNA or RNA and only them
        assert type == "DNA" or "RNA", type
        self.type = type
        #check if error type smaller than 1
        assert error_rate <= 1, error_rate
        self.error_rate = error_rate
        self.flag = 0

    # This func get a string a return the reversed string
    def reverse(self, string):
        rev_str = ""
        #currently char
        for char in string:
            # Adding the current char in front of the rest string
            rev_str = char + rev_str
        return rev_str

    # this func return the complementary base for a base
    def couples_dict(self, char):
        #if its dna polymerase
        if self.type == "DNA":
            dictionary = {'a': 'U', 'A': 'U', 'c': 'G', 'C': 'G', 'g': 'C', 'G': 'C', 't': 'A', 'T': 'A',}
        # if its rna polymerase
        else:
            dictionary = {'a': 'T', 'A': 'T', 'c': 'G', 'C': 'G', 'g': 'C', 'G': 'C', 'u': 'A', 'U': 'A',}
        return dictionary[char]

    # This func get a sequence as list and a place for a mutation and choose randomly the mutated base
    def mutation(self, seq, place):
        # bases options for dna seq
        if self.type == "DNA":
            mutation_bases = ['A', 'C', 'G', 'T']
        # bases options for dna seq
        elif self.type == "RNA":
            mutation_bases = ['A', 'C', 'G', 'U']
        # removing the current base
        mutation_bases.remove(seq[place])
        # choosing the mutated base
        mut_base = random.choice(mutation_bases)
        return mut_base

    # This func get a dna sequence and return it's complementary rna sequence as a string
    # the flag is for knowing if should return complimentary seq or same as func get (for creating mutation only)
    def transcribe(self, dna_seq):
        dna_upper_seq = dna_seq.upper()
        current_seq = []
        current_seq[:0] = dna_upper_seq
        # first, check if the sequence match the type of polymerase (dna in this case)
        if self.type == "DNA":
            assert check_dna(dna_seq), dna_upper_seq
        # first, check if the sequence match the type of polymerase (rna in this case)
        elif self.type == "RNA":
            assert check_rna(dna_seq), dna_upper_seq
        seq_length = len(dna_upper_seq)
        mut_num = math.ceil(float(seq_length * self.error_rate))
        # choose random places for the mutations
        mut_places = random.sample(range(seq_length), mut_num)
        # creating the mutation in every place in mut_places
        for place in mut_places:
            current_seq[place] = self.mutation(current_seq, place)
        str_seq= ''.join(current_seq)
            # continue to complementary seq
        if self.flag == 0:
            new_seq = ""
            # Adding rna bases to the new seq accordingly to their complementary dna bases
            for i in str_seq:
                new_seq = new_seq + self.couples_dict(i)
            # Changing the order for 5' to 3'
            rna_seq = self.reverse(new_seq)
            return rna_seq
        # return the duplicate seq but the mutations created
        elif self.flag == 1:
            return str_seq

#this class represent a ribosome
class Ribosome:
    #initial a ribosome
    def __init__(self, genetic_code, start_codons):
        self.genetic_code = genetic_code
        self.start_codons = start_codons

    # This func get a sequence and cut away the remainder of division by 3, and return it as a string
    def end_by_triplets(self, sequence):
        # Calculate the remainder
        remainder = len(sequence) % 3
        # Set the end according the triplets
        end = len(sequence) - remainder
        return end

    # This func get a sequence and create a reading frame according to start codon and stop codons if there are any
    # The func will return the longest reading frame possible as a string
    def framing(self, sequence):
        stop_codons = []
        # getting the stop codons from the genetic code
        for key in self.genetic_code.keys():
            if self.genetic_code[key] == None:
                stop_codons.append(key)
        # Final frame will be here
        frame = ""
        # Current frame we check
        temp_seq = ""
        # Checking if there is a start codon
        for i in range(0,len(sequence)-2):
            sub_start = sequence[i:i+3]
            # if this id a start codon
            if sub_start in self.start_codons:
                start = i
                # Setting the start of the temp seq according to the start codon
                temp_seq = sequence[start:len(sequence)]
                # Setting the end of the frame according to the triplets
                end = self.end_by_triplets(temp_seq)
                # Setting the end of the frame according to the stop codons if any
                for j in range(3, len(temp_seq) + 1, 3):
                    # sub - current triple
                    sub = temp_seq[j:j + 3]
                    if sub in stop_codons:
                        end = j
                # Cutting the temp sequence from start codon (above) to the final end
                temp_seq = temp_seq[0:end]
                # Setting the longest frame possible
                if len(temp_seq) > len(frame):
                    frame = temp_seq
        # if there is no reading frame return None
        if frame == "":
            return None
        return frame

    # This func get a rna sequence and translate it to codons of amino acid
    # The func will return a list of the codons
    def translate(self, rna_seq):
        # Set the reading frame of the rna sequence that the function got
        reading_frame = self.framing(rna_seq)
        codons_list = []
        #if there is no reading frame (and no codons)
        if reading_frame is None:
            return None
        else:
            #adding codons to the list
            for i in range(0,len(reading_frame)-2, 3):
                codon = reading_frame[i:i+3]
                codons_list.append(codon)
        return codons_list

    #this function get a rna sequence and translate it, return the protein sequence
    def synthesize(self, rna_seq):
        #first, check if the rna seq is actually a ren seq
        assert check_rna(rna_seq), rna_seq
        protein = ""
        temp = ""
        #translating the codons to amino acids
        codons_list = self.translate(rna_seq)
        #adding every amino acid to the protein seq
        if codons_list is None:
            protein = "NA"
            return protein
        for codon in codons_list:
            temp = self.genetic_code[codon]
            if temp:
                protein = protein + temp
        return protein

#this class represent a cell
class Cell:
    #initial a cell
    def __init__(self, name, genome, num_copies, genetic_code, start_codons, division_rate):
        self.name = name
        #check if the genome actually contains correct dna sequence as genes
        for gene in genome:
            assert check_dna(gene), gene
        self.genome = genome
        # check if num of copies of the genome is actually bigger than 0 and integer
        assert num_copies > 0, num_copies
        if type(num_copies) != int:
            assert num_copies.is_integer() == True, num_copies
        self.num_copies = num_copies
        self.genetic_code = genetic_code
        self.start_codons = start_codons
        # check if division rate is actually bigger than 1 and integer
        assert division_rate > 1, division_rate
        if type(division_rate) != int:
            assert division_rate.is_integer() == True, division_rate
        self.division_rate = division_rate
        self.dna_pol = Polymerase("DNA")
        self.rna_pol = Polymerase("RNA")
        self.rib = Ribosome(genetic_code, start_codons)

    #printing the cell info
    def print_cell(self):
        print("Original cell: " + "<{}, {}, {}>".format(self.name, self.num_copies, self.division_rate))

    # This func get a sequence, index of the triplets repeats start, and the size of the suc string (between 1-6)
    # The func check if there are more than 3 repeats and at the end return the number of them
    def how_many(self, dna_seq, index, size):
        # At this point, we already got a triple
        counter = 3
        # Setting the next place we want to check
        next_first_place = index + (3 * size)
        # Checking how much ths sub-string repeats till the end of sequence
        for x in range(next_first_place, len(dna_seq) - size + 1, size):
            y = x + size
            # If there is a match, increase the counter by 1
            if dna_seq[x:y] == dna_seq[index:index + size]:
                counter = counter + 1
            else:
                return counter
        return counter

    # This func get a dna sequence, and check if there are simple (1-6) strings repeatedly 3+ times at the sequence
    # If there are, the func return a dict with the repeats and the max num of appears in the seq
    def find_srr(self, dna_seq):
        #first, check if the dna seq is actually a dna seq
        assert check_dna(dna_seq), dna_seq
        repeats = {}
        temp_sub = ""
        # Size of the sub-string we check (i is for the left limit, j is for the right limit)
        for size in range(1, 7):
            i = 0
            # As long the sub string didnt cross the end of the sequence
            while i < len(dna_seq) - size + 1:
                j = i + size
                # If the current sub-string is similar to the temp sub-string we already  check, keep moving
                if dna_seq[i:j] == temp_sub or dna_seq[i:j] in repeats.keys():
                    i = i + 1
                # If not, check if there is triple of the same sub-string
                elif dna_seq[i:j] == dna_seq[i + size:j + size] and dna_seq[i:j] == dna_seq[i + (2 * size):j + (2 * size)]:
                    # Updating the temp sub-string accordingly to the triple we found
                    temp_sub = dna_seq[i:j]
                    # Calling a func to find out how many times the sub-string appears
                    counter = self.how_many(dna_seq, i, size)
                    # Updating the repeats string accordingly
                    repeats[temp_sub] = counter
                    # Moving on
                    i = i + 1
                # If there isn't any repeated sun string now, reset the temp sub string and moving on the next check
                else:
                    temp_sub = ""
                    i = i + 1
        repeats = {key: val for key, val in sorted(repeats.items())}
        #if there is no repeats
        if repeats == {}:
            return None
        return repeats

    # create a list of division cells
    def mitosis(self):
        #duplicate the cell according to division rate of the cell
        cell_list = [self] * self.division_rate
        return cell_list

    # for the second cell from meiosis return the complementary sequence, replacing the base and reversing the seq
    def complementary(self, dna_seq):
        com_seq = ""
        dictionary = {'a': 'T', 'A': 'T', 'c': 'G', 'C': 'G', 'g': 'C', 'G': 'C', 't': 'A', 'T': 'A', }
        # for every base, chang to the complementary base
        for char in dna_seq:
            com_seq = dictionary[char] + com_seq
        return com_seq

    # this func devide a cell into two cells, one with the original genome and one with the complementary genome
    def meiosis(self):
        # if cant do meiosis (odd num of copies)
        if self.num_copies % 2 != 0 :
            return None
        # create 2 cells, one with the original genome and one with the complementary genome
        else:
            first = Cell(self.name, self.genome, self.num_copies / 2, self.genetic_code, self.start_codons, self.division_rate)
            com_genome = []
            # create the complementary genome
            for gene in self.genome:
                com_gene = self.complementary(gene)
                com_genome.append(com_gene)
            second = Cell(self.name, com_genome, self.num_copies / 2, self.genetic_code, self.start_codons, self.division_rate)
            cell_list = [first, second]
            return cell_list

    #this func return the repertoire of every dna sequence in the genome, contains the srr, rna sequence, and protein sequence
    def repertoire(self):
        rep_list = []
        # for every seq from the genome
        for gene in self.genome:
            # getting the rna sequence
            self.dna_pol.flag = 0
            rna = self.dna_pol.transcribe(gene)
            # getting the protein sequence if there is one coded by the rna
            if self.rib.translate(rna) == None:
                coding = "Non-coding RNA"
            else:
                coding = self.rib.synthesize(rna)
            print_srr = ""
            # getting the srr found in the dna sequence if there are
            if self.find_srr(gene) == None:
                srr = "No simple repeats in DNA sequence"
                t = (srr, rna, coding)
            else:
                srr = self.find_srr(gene)
                # getting the srr into a string
                for key in srr.keys():
                    print_srr = print_srr + key + "," + str(srr[key]) + ";"
                print_srr = print_srr[0:-1]
                # make a touple for every gene (dna seq)
                t = (print_srr, rna, coding)
            rep_list.append(t)
        return rep_list

#this class represent a prokaryotic cell
class ProkaryoticCell(Cell):
    #initial a prokaryotic cell
    def __init__(self, name, genome):
        self.name = name
        # check if the genome actually contains correct dna sequence as genes
        for gene in genome:
            assert check_dna(gene), gene
        self.genome = genome
        self.num_copies = 1
        self.genetic_code = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA': None, 'UAG': None,
        'UGC':'C', 'UGU':'C', 'UGA': 'U', 'UGG':'W'}
        self.start_codons =  ["AUG", "GUG", "UUG"]
        self.division_rate = 4
        self.dna_pol = Polymerase("DNA")
        self.rna_pol = Polymerase("RNA")
        self.rib = Ribosome(self.genetic_code, self.start_codons)

#this class represent a eukaryotic cell
class EukaryoticCell(Cell):
    # initial a eukaryotic cell
    def __init__(self, name, genome):
        self.name = name
        # check if the genome actually contains correct dna sequence as genes
        for gene in genome:
            assert check_dna(gene), gene
        self.genome = genome
        self.num_copies = 2
        self.genetic_code = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA': None, 'UAG': None,
        'UGC':'C', 'UGU':'C', 'UGA': None, 'UGG':'W'}
        self.start_codons = ["AUG"]
        self.dna_pol = Polymerase("DNA")
        self.rna_pol = Polymerase("RNA")
        self.rib = Ribosome(self.genetic_code, self.start_codons)

#this class represent a neuron cell
class NeuronCell(EukaryoticCell):
    # initial a neuron cell
    def __init__(self, name, genome):
        #inheritance of initialization
        EukaryoticCell.__init__(self, name, genome)
        self.division_rate = 2

#this class represent a stem cell
class StemCell(EukaryoticCell):
    # initial a stem cell
    def __init__(self, name, genome):
        # inheritance of initialization
        EukaryoticCell.__init__(self, name, genome)
        self.division_rate = 3

#this class represent a mutant cell
class MutantCell(StemCell):
    #initial a mutant cell
    def __init__(self, genome, num_mutations=0, error_rate=0.05):
        StemCell.__init__(self, "MutantCell", genome)
        self.name = "MutantCell"
        self.dna_pol.error_rate = error_rate
        self.mutations = num_mutations

    # this function create a list of division cells, if the cell have more than 10 mutation cell become a cancer cell
    def mitosis(self):
        cells_list = [self]
        mut_genome = []
        mutation_num = self.mutations
        #create the mutated genome for daughter cells
        for gene in self.genome:
            self.dna_pol.flag = 1
            mut_gene = self.dna_pol.transcribe(gene)
            mut_genome.append(mut_gene)
            mutation_num += math.ceil(float(len(gene) * self.dna_pol.error_rate))
        # for more than 10 mutations, the cell become a cancer cell
        if mutation_num > 10:
            cancer_cell = CancerCell(mut_genome, mutation_num)
            # duplicate according to the division rate of the cell
            for dev in range(1, self.division_rate):
                cells_list.append(cancer_cell)
        # duplicate according to the division rate of the cell
        else:
            mut_cell = MutantCell(mut_genome, mutation_num)
            for dev in range(1, self.division_rate):
                cells_list.append(mut_cell)
        return cells_list

#this class represent a cancer cell
class CancerCell(MutantCell):
    #initial a mutant cell
    def __init__(self, genome, num_mutations=0):
        MutantCell.__init__(self,genome, num_mutations)
        self.name = "CancerCell"
        self.division_rate = 10

# this class represent a calssifier for pattern, proteins and domains
class SequenceClassifier:
    def __init__(self, pattern_file):
        self.pattern_dict = self.__patterns_to_domains(pattern_file)

    # this function translate a pattern to regex
    def pattern_to_regex(self, pattern):
        reg = ""
        replace_dict = {'[' : '[', ']' : ']', '-' : '', 'x' : '.', '(' : '{', ')' : '}',
                        '{' : '[^', '}' : ']', '<' : '\A', '>' : '\Z'}
        # replacing each char to the matching regex char
        for char in pattern:
            if char in replace_dict:
                reg += replace_dict[char]
            elif 65 <= ord(char) <= 90 or 48 <= ord(char) <= 57:
                reg += char
            elif 97 <= ord(char) <= 122:
                reg += chr(ord(char) - 32)
           # if there is chat that doesnt match regex expression, return empty string
            else:
                return ""
        return reg

    #this function get a dict of patterns to domains and create a dict of regex expressions and domains
    def __prosite_to_python(self, pattern_dict):
        reg_dict = {}
        temp_reg = ""
        # create a regex expression for every pattern
        for pattern in pattern_dict:
            temp_reg = self.pattern_to_regex(pattern)
            # check if the expression is valid regex
            try:
                re.compile(temp_reg)
            except:
                print("ValueError" + temp_reg)
            else:
                reg_dict[temp_reg] = pattern_dict[pattern]
        return reg_dict

    # this function get a csv file with patterns and their domains, and return a dict of regex expression and their domains
    def __patterns_to_domains(self, pattern_file):
        pattern_dict = {}
        reg_dict = {}
        # open the file and copy the information to a dict
        assert os.path.exists(pattern_file), pattern_file
        with open(pattern_file, 'r') as f:
            data = csv.reader(f)
            for row in data:
                temp_key = row[0]
                temp_val = row[1]
                pattern_dict[temp_key] = temp_val
        # replacing the dict from patterns to regex by previous function
        reg_dict = self.__prosite_to_python(pattern_dict)
        return reg_dict

    # this function get a gene list and return a protein list
    def gene_to_protein(self, gene_list):
        seq_list = []
        # createin an eukaryotic gene for using his polymerase anf ribosome to convert every gene to protein
        current_cell = EukaryoticCell("EukaryoticCell", gene_list)
        # for every seq from the genome
        for gene in current_cell.genome:
            # getting the rna sequence
            current_cell.dna_pol.flag = 0
            rna_seq = current_cell.dna_pol.transcribe(gene)
            # getting the protein sequence if there is one coded by the rna
            protein_seq = current_cell.rib.synthesize(rna_seq)
            if protein_seq is not None:
                seq_list.append(protein_seq)
        return seq_list

    # this function get a list of genomic sequences and return dict for each protein from the gen seq and his domains
    def find_domains(self, seq_list):
        # call a func that convert the gen seq list to proteins list
        protein_seq_list = self.gene_to_protein(seq_list)
        temp_dict = {}
        temp_domains = []
        # check if protein is valid
        for seq in protein_seq_list:
            if seq == "NA":
                continue
            # if there is a match between the protein seq and this class patterns dict, finding the domains
            for pattern in self.pattern_dict:
                m = re.match(pattern, seq)
                if m:
                    domain = self.pattern_dict[pattern]
                    # assuring that every domain appears once
                    if domain not in temp_domains:
                        temp_domains.append(self.pattern_dict[pattern])
            # if there is no domain, value will be NA - Not Available
            if len(temp_domains) == 0:
                temp_dict[seq] = ["NA"]
            # adding the domains
            else:
                domains = [';'.join(temp_domains)]
                temp_dict[seq] = domains
            temp_domains = []
        return temp_dict

    # this function get a list of genomic sequences and path for csv file, and create clasified dict for protein seq to domain
    def classify(self, seq_list, csv_file):
        # calling for a func to convert each gen seq to protein seq and his domain
        classifies_dict = self.find_domains(seq_list)
        # create the csv file
        assert os.path.exists(csv_file), csv_file
        with open(csv_file, 'w') as f:
            csvwriter = csv.writer(f)
            csvwriter.writerow(['Sequence', 'Domains'])
            for seq in classifies_dict:
                csvwriter.writerow([seq, classifies_dict[seq]])

# this function make a simulation for specific values for 3 repeats and return the results in a list
def simulation(genome, dev_num, error_rate):
    info = []
    cancer_cells_num = 0
    mutant_cells_num = 0
    diff_proteins_num = 0
    # 3 repeats
    for i in range(3):
        temp_dev = dev_num
        temp_genome = []
        current_cell = MutantCell([str(genome)], 0, error_rate)
        cells_list = [current_cell]
        new_list = []
        index = 0
        # while there are more cycles left
        while temp_dev > 0:
            types = []
            for cell in cells_list:
                if cell not in types:
                    types.append(cell)
            #for every cell in the cell list making mitosis according to the cell type (same type doing the same mitosis)
            for celltype in types:
                for cell in cells_list:
                    if cell == celltype:
                        index += 1
                for num in range(index):
                    new_list.extend(celltype.mitosis())
                index = 0
            cells_list = new_list
            new_list = []
            temp_dev = temp_dev - 1
        # in the final cell list, counting the number of mutant and cancer cells
        for cell in cells_list:
            temp_genome.extend(cell.genome)
            if cell.name == "MutantCell":
                mutant_cells_num = mutant_cells_num + 1
            elif cell.name == "CancerCell":
                cancer_cells_num = cancer_cells_num + 1
        temp_diff_genome = list(set(temp_genome))
        diff_proteins = set()
        c = EukaryoticCell("c", temp_diff_genome)
        # counting the number of unique proteins in the genome
        for info in c.repertoire():
            if info[-1] != "Non-coding RNA":
                diff_proteins.add(info[-1])
        diff_proteins_num = diff_proteins_num + len(diff_proteins)
    # calculating the average of the 3 repeats results
    mutant_cells_num = mutant_cells_num / 3
    cancer_cells_num = cancer_cells_num / 3
    diff_proteins_num = diff_proteins_num / 3
    info = [mutant_cells_num, cancer_cells_num, diff_proteins_num]
    return info


random.seed(1)
assert len(sys.argv) > 1, "did not get file"
fasta_file = sys.argv[1]
assert fasta_file.endswith(".fa"), fasta_file
assert os.path.exists(fasta_file), fasta_file
genomes = []
with open(fasta_file) as file:
    for seq_record in seqIO.parse(file, "fasta"):
        genomes.append(seq_record)
error_rates = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
with open("exercise4_318731411_316393974.csv", 'w') as file:
    csvwriter = csv.writer(file)
    csvwriter.writerow(['Sequence Name','Error Rate','Number of Cycles','Mutant Cells','Cancer Cells','Unique Proteins'])
    for genome in genomes:
        cancer_cells_y = []
        mutated_cells_y = []
        protein_seq_num_y_5 = []
        protein_num_from_all_dev = []
        array_of_variables = []
        for dev_num in range(1, 6):
            for error_rate in error_rates:
                sim_info = simulation(genome.seq, dev_num, error_rate)
                if dev_num == 5:  # taking some data from the 5th division
                    cancer_cells_y.append(sim_info[1])  # taking the data for 5A
                    mutated_cells_y.append(sim_info[0])  # taking the data for 5A
                    protein_seq_num_y_5.append(sim_info[2])  # taking the data for 5A
                protein_num_from_all_dev.append(sim_info[2])  # taking the data for 5B1 and 5B2
                csvwriter.writerow([genome.name, error_rate, dev_num, sim_info[0], sim_info[1], sim_info[2]])
        array_of_variables.append(cancer_cells_y)  # next 3 lines create ordered array of Y in graph 5A
        array_of_variables.append(mutated_cells_y)
        array_of_variables.append(protein_seq_num_y_5)
        counter = 1
        names_of_y = {1: 'cancer cells', 2: 'mutant cells', 3: 'num of proteins'}  # dict of Y labels
        for array in array_of_variables:  # creating a meta graph of all three graphs of parameters per error rate
            graph_x = error_rates
            graph_y = array
            plt.subplot(3, 1, counter)
            plt.plot(graph_x, graph_y)
            y_string = names_of_y[counter]
            plt.title(y_string)
            plt.xlabel("Error rate")
            plt.ylabel(y_string)
            counter += 1
        plt.suptitle(genome.name)  # naming the meta graph with the name of the seq, so we know which gene it is
        plt.savefig("exercise4_318731411_316393974_5A_{}.png".format(genome.name))  # saving it in a file
        plt.show()  # showing it
        plt.clf()  # for the terminal
        # 5b1:
        plt.plot(error_rates, protein_num_from_all_dev[0:11], 'r--',  # creating 5 plots of proteins per error rate
                 error_rates, protein_num_from_all_dev[11:22], 'bs',  # for each division
                 error_rates, protein_num_from_all_dev[22:33], 'g^',
                 error_rates, protein_num_from_all_dev[33:44], 'k',
                 error_rates, protein_num_from_all_dev[44:55], 'm:')
        plt.title(genome.name)  # name of the gene
        plt.ylabel("num of proteins")
        plt.xlabel("error rate")
        plt.legend(["one division", "two divisions", "three divisions", "four divisions",
                    "five divisions"], loc="upper left")  # creating a legend that would explain each plot
        plt.savefig("exercise4_318731411_316393974_5B1_{}.png".format(genome.name))  # saving it in file
        plt.show()  # showing it
        plt.clf()  # for the terminal
        # 5b2: a plot of sum of diff proteins of all error rates per num of cycle
        num_of_proteins_1 = 0
        for num in range(11):
            num_of_proteins_1 += protein_num_from_all_dev[num]  # sum of protein per division number
        num_of_proteins_2 = 0
        for num in range(11):
            num_of_proteins_2 += protein_num_from_all_dev[num + 11]
        num_of_proteins_3 = 0
        for num in range(11):
            num_of_proteins_3 += protein_num_from_all_dev[num + 22]
        num_of_proteins_4 = 0
        for num in range(11):
            num_of_proteins_4 += protein_num_from_all_dev[num + 33]
        num_of_proteins_5 = 0
        for num in range(11):
            num_of_proteins_5 += protein_num_from_all_dev[num + 44]
        x = [1, 2, 3, 4, 5]  # list of divisions
        y = [num_of_proteins_1, num_of_proteins_2, num_of_proteins_3, num_of_proteins_4, num_of_proteins_5]
        plt.scatter(x, y)  # creating a scatter plot
        plt.title(genome.name)  # the name of the gene which has this data
        plt.ylabel("num of proteins")
        plt.xlabel("divisions")
        plt.savefig("exercise4_318731411_316393974_5B2_{}.png".format(genome.name))  # saving it in a file
        plt.show()  # showing the graph
        plt.clf()  # for the terminal
