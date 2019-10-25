#!/usr/bin/python
# -*- coding: UTF-8 -*-

import argparse
import sys
import re

# Usage & optional inputs
frame_choices = ['+1', '+2', '+3', '-1', '-2', '-3']
parser = argparse.ArgumentParser(description='An ORF Finder, genome sequence files and output file are in the same folder as this programme. Every input fasta file could contain one or more nucleatide sequences. Enjoy! :)')
parser.add_argument('-i', dest='input_files', default='input.fasta', help='input filenames, use \',\' to separate multiple files')
parser.add_argument('-o', dest='output_file', default='proteins.fasta', help='output filename, only one')
parser.add_argument('-f', dest='frame_number', help='restrict to one single frame', choices=frame_choices)
parser.add_argument('-s', dest='min_size', default=50, type=int, help='define the mininum amino acids for each ORF')
args = parser.parse_args()
read_files = args.input_files.split(',')
write_file = args.output_file
min_size = int(args.min_size)
frame_number = args.frame_number

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)

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

# Print optional information
def printInputInfo():
    input_file_number = len(read_files)
    if input_file_number == 1:
        print('Input file:', read_files[0])
    else:
        print('Input files', ' '.join(read_files))
    print('Output file:', write_file)
    print('Minimun Amino Acids:', min_size)
    if frame_number:
        print('Restrict to frame: %s' % frame_number)
    else:
        print('No valid frame number specified, will find for all 6 frames')

# Fasta file check
def fastaFormatCheck(file):
    content = list(filter(lambda l: l != '\n', file.readlines()))
    if not content[0].startswith('>') or content[-1].startswith('>'):
        raise TypeError('This is not a valid fasta file.')
    idx = 0
    trig = False
    while idx < len(content):
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

# Read all input files, check the validity of them and format all sequences in 1 list
def readAllSeqs(read_files):
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
        except TypeError as errObj:
            print('{0}: {1}'.format(file, errObj))
            sys.exit(1)
        except:
            print('Cannot find file \'%s\'. Please input exist filename.' % file)
            sys.exit(2)
        finally:
            f.close()
    return all_contents

# Reverse Complement and List All or Single One frames
def reverseComplement(seq):
    transTab = str.maketrans('ATGC', 'TACG')
    comp = seq.translate(transTab)
    revcamp = comp[::-1]
    return revcamp

def fillAll6(seq):
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
    idx = int(frame_number[1]) - 1
    if '+' in frame_number:
        return {frame_number: seq[idx:]}
    else:
        rc_version = reverseComplement(seq)
        return {frame_number: rc_version[idx:]}

# translation
def proteinTranslation(seq):
    amino_acids = ''
    seq_list = [seq[i:i+chunk_size] for i in range(0, len(seq), chunk_size)]
    for dna in seq_list:
        if len(dna) < 3:
            break
        if 'N' in dna:
            amino_acids += 'X'
            continue
        amino_acids += DNA_codon[dna]
    return amino_acids

def translateFromFrames(frames):
    translated = {}
    for key in frames:
        trans = proteinTranslation(frames[key])
        translated[key] = trans
    return translated

# find the longest ORF
def findLongestORF(seq, frame = '+1'):
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

# main entrance
def find(seq):
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
    output = open(write_file, 'w')
    all_seqs = readAllSeqs(read_files)
    for item in all_seqs:
        find(item)
    output.close()