#!/usr/bin/python
# -*- coding: UTF-8 -*-

import argparse
import sys
import re

# Usage & optional inputs
frame_choices = ['+1', '+2', '+3', '-1', '-2', '-3']
parser = argparse.ArgumentParser(description='An ORF Finder, genome sequence files and output file are in the same folder as this programme. Every input fasta file could contain one or more nucleotide sequences. It will skip all gaps and remove all ORFs with unknown codons automatically.')
parser.add_argument('-i', dest='input_files', required=True, help='input filenames, use \',\' to separate multiple files')
parser.add_argument('-o', dest='output_file', default='proteins.fasta', help='output filename, the default output file is proteins.fasta')
parser.add_argument('-f', dest='frame_number', help='restrict to one single frame, if not specified, will produce for all frames', choices=frame_choices)
parser.add_argument('-s', dest='min_size', default=50, type=int, help='define the minimum amino acids for each ORF, the default size is 50')

if 'man' in sys.argv: # python3 group2_task1.py man, simulating the man page in Linux, make sure the man page is under the same folder with this script and run this script under its own folder
        manpage = open('group2_task1.txt', 'r')
        print(''.join(manpage.readlines()))
        manpage.close()
        sys.exit(0)
try:
    if '-h' in sys.argv:
        sys.exit(0)
    args = parser.parse_args()
except:
    parser.print_help() # anyway, the programme will output the whole help information :)
    sys.exit(1)
read_files = args.input_files.split(',')
write_file = args.output_file
min_size = int(args.min_size)
frame_number = args.frame_number

# Global Variables
DNA_codon = {
        'AAA' : 'K', 'AAC' : 'N', 'AAG' : 'K', 'AAT' : 'N',
        'ACA' : 'T', 'ACC' : 'T', 'ACG' : 'T', 'ACT' : 'T',
        'AGA' : 'R', 'AGC' : 'S', 'AGG' : 'R', 'AGT' : 'S',
        'ATA' : 'I', 'ATC' : 'I', 'ATG' : 'M', 'ATT' : 'I',
        'CAA' : 'Q', 'CAC' : 'H', 'CAG' : 'Q', 'CAT' : 'H',
        'CCA' : 'P', 'CCC' : 'P', 'CCG' : 'P', 'CCT' : 'P',
        'CGA' : 'R', 'CGC' : 'R', 'CGG' : 'R', 'CGT' : 'R',	
        'CTA' : 'L', 'CTC' : 'L', 'CTG' : 'L', 'CTT' : 'L',
        'GAA' : 'E', 'GAC' : 'D', 'GAG' : 'E', 'GAT' : 'D',
        'GCA' : 'A', 'GCC' : 'A', 'GCG' : 'A', 'GCT' : 'A',
        'GGA' : 'G', 'GGC' : 'G', 'GGG' : 'G', 'GGT' : 'G',
        'GTA' : 'V', 'GTC' : 'V', 'GTG' : 'V', 'GTT' : 'V',
        'TAA' : '*', 'TAC' : 'Y', 'TAG' : '*', 'TAT' : 'Y',
        'TCA' : 'S', 'TCC' : 'S', 'TCG' : 'S', 'TCT' : 'S',
        'TGA' : '*', 'TGC' : 'C', 'TGG' : 'W', 'TGT' : 'C',
        'TTA' : 'L', 'TTC' : 'F', 'TTG' : 'L', 'TTT' : 'F'  
    }
chunk_size = 3

def printInputInfo():
    '''
        Print all optional information
    '''
    input_file_number = len(read_files)
    if input_file_number == 1:
        print('Input file:', read_files[0])
    else:
        print('Input files', ' '.join(read_files))
    print('Output file:', write_file)
    print('Minimum Amino Acids:', min_size)
    if frame_number:
        print('Restrict to frame: %s' % frame_number)
    else:
        print('No valid frame number specified, will find for all 6 frames')

def fastaFormatCheck(file):
    '''
        Fasta file format check
    '''
    content = list(filter(lambda l: l != '\n', file.readlines())) # skip all null lines
    if not content[0].startswith('>') or content[-1].startswith('>'):
        raise TypeError('This is not a valid fasta file.')
    idx = 0
    trig = False
    while idx < len(content):
        if not content[idx].startswith('>') and re.search(r' ', content[idx]):
            print('Gap is found at line {0} in {1}. Will remove it automatically.'.format(idx + 1, file.name))
            content[idx] = re.sub('[\s+]', '', content[idx])
        if trig:
            if content[idx].startswith('>'):
                raise TypeError('This is not a valid fasta file.')
            else:
                trig = False
        else:
            if content[idx].startswith('>'):
                trig = True
        idx += 1
    return content

def readAllSeqs(read_files):
    '''
        Read all input files, check the validity of them and format all sequences in 1 list
    '''
    all_contents = []
    for file in read_files:
        try:
            f = open(file, 'r')
            lines = fastaFormatCheck(f)
            for line in lines:
                if line.startswith('>'):
                    all_contents.append(line.split()[0].replace('>', '') + ':')
                else:
                    all_contents[-1] += line.rstrip('\n').upper()
            f.close()
        except TypeError as errObj:
            print('{0}: {1}'.format(file, errObj))
            sys.exit(2)
        except:
            print('Cannot find file \'%s\'. Please input exist filename.' % file)
            sys.exit(3)
    return all_contents

def reverseComplement(seq):
    '''
        Reverse Complement and List All or Single One frames
    '''
    transTab = str.maketrans('ATGC', 'TACG')
    comp = seq.translate(transTab)
    revcamp = comp[::-1]
    return revcamp

def fillAll6(seq):
    '''
        for all 6 frames
    '''
    all_frames = {}
    all_frames['+1'] = seq
    all_frames['+2'] = seq[1:]
    all_frames['+3'] = seq[2:]
    rc_version = reverseComplement(seq)
    all_frames['-1'] = rc_version
    all_frames['-2'] = rc_version[1:]
    all_frames['-3'] = rc_version[2:]
    return all_frames

def fillOneFrame(seq, frame_number):
    '''
        for 1 single frame
    '''
    idx = int(frame_number[1]) - 1
    if '+' in frame_number:
        return {frame_number: seq[idx:]}
    else:
        rc_version = reverseComplement(seq)
        return {frame_number: rc_version[idx:]}

def proteinTranslation(seq):
    '''
        translation
    '''
    amino_acids = ''
    seq_list = [seq[i:i+chunk_size] for i in range(0, len(seq), chunk_size)]
    for dna in seq_list:
        if len(dna) < 3:
            break
        if re.search(r'[^ACGT]', dna):
            amino_acids += 'X'
            continue
        amino_acids += DNA_codon[dna]
    return amino_acids

def translateFromFrames(frames):
    '''
        translate into proteins for all sequences in the frames
    '''
    translated = {}
    for key in frames:
        trans = proteinTranslation(frames[key])
        translated[key] = trans
    return translated

def findLongestORF(seq, frame = '+1'):
    '''
        find the longest ORF in 1 sequence
    '''
    ORFs = {}
    tmp = ''
    idx = 0
    if '2' in frame:
        idx = 1
    elif '3' in frame:
        idx = 2
    for char in seq:
        if tmp:
            if char != '*':
                tmp += char
            else:
                if (len(tmp) >= min_size) & ('X' not in tmp):
                    ORFs[idx - len(tmp)* 3 + 1] = tmp
                tmp = ''
        elif char == 'M':
            tmp += char
        idx += 3
    if (len(tmp) >= min_size) & ('X' not in tmp):
        ORFs[idx - len(tmp)* 3 + 1] = tmp
    return ORFs

def find(seq):
    '''
        main function
    '''
    name_and_seq = seq.split(':')
    name_start = name_and_seq[0]
    seq = name_and_seq[1]
    all_length = len(seq)
    frames = {}
    if (frame_number):
        frames = fillOneFrame(seq, frame_number)
    else:
        frames = fillAll6(seq)
    translated_frames = translateFromFrames(frames)
    for key in translated_frames:
        ORFs = findLongestORF(translated_frames[key], key)
        if len(ORFs):
            index = 1
            for start, orf in ORFs.items():
                if '-' in key:
                    start = all_length - start + 1
                output.write('>%s_F%s_%04d\t%s\t%d\t%d\n' % (name_start, key[1], index, key, len(orf), start))
                output.write('\n'.join([orf[i:i+60] for i in range(0,len(orf),60)]) + '\n')
                index += 1
    print('ORFs have been found for: %s, please check in output file: %s' % (name_start, write_file))

if __name__ == '__main__':
    printInputInfo()
    all_seqs = readAllSeqs(read_files)
    output = open(write_file, 'w')
    for item in all_seqs:
        find(item)
    output.close()