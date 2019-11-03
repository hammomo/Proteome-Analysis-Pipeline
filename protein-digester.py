#!/usr/bin/python
# -*- coding: UTF-8 -*-

import argparse
import sys
import re

# Usage & optional inputs
enzymes = ['t', 'k', 'r', 'e']
parser = argparse.ArgumentParser(description='This is a protein digester, using the user-defined choice to digest proteins into small peptide sequences. Enzyme choice could be selected from "t. Trypsin, k. Endoproteinase Lys-C, r. Endoproteinase Arg-C, e. V8 proteinase (Glu-C)". Please input the number of these choices.')
parser.add_argument('-e', '--enzyme', help='select one digesting enzyme to cut proteins', choices=enzymes)
parser.add_argument('-i', '--input', default='proteins.fasta', help='input file which contains protein sequences')
parser.add_argument('-o', '--output', default='peptides.peps', help='output file which contains the result of this programme (peptides)')
parser.add_argument('-c', '--cleavage', type=int, default=0, help='input the number of missed cleavages')
args = parser.parse_args()
enzyme = args.enzyme
input_file = args.input
output_file = args.output
cleavage_number = int(args.cleavage)

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)

if not enzyme:
    print('Enzyme shoule be specified!')
    parser.print_help()
    sys.exit(3)

if not input_file:
    print('Input file must be specified!\n')
    parser.print_help()
    sys.exit(1)

enzyme_rules = {
    # 't': r'([RK])([^P])',
    # 'k': r'(K)([^P])',
    # 'r': r'(R)([^P])',
    # 'e': r'(E)([^P])'
    't': r'([RK])',
    'k': r'(K)',
    'r': r'(R)',
    'e': r'(E)'
}

enzyme_dict = {
    't': 'trypsin',
    'k': 'lys-c',
    'r': 'arg-c',
    'e': 'glu-c'
}

def printInputInfo():
    print('Input filename:', input_file)
    print('Output filename:', output_file)
    print('Enzyme digesting choice:', enzyme_dict[enzyme])
    if cleavage_number:
        print('The number of missed cleavages:', cleavage_number)
    else:
        print('No missed cleavage')

def enzymeFullCut(seq, rule):
    splits = re.split(rule, seq)
    result = [''.join([splits[i], splits[i+1]]) for i in range(0, len(splits) - 1, 2)]
    if splits[-1]:
        result.append(splits[-1])
    i = 0
    while i < len(result):
        if result[i].startswith('P') and i != 0:
            result[i-1] += result[i]
            result.pop(i)
        else:
            i += 1
    # result = re.sub(rule, r'\1 \2', seq).split()
    return result

def missCleavage(peps, number):
    i = 0
    result = []
    while i < len(peps) - number:
        tmp = ''
        for j in range(0, number + 1):
            tmp += peps[i+j]
        result.append(tmp)
        i += 1
    return result

def readFile(file):
    proteins = {}
    try:
        fileObj = open(file, 'r')
        name = ''
        for line in fileObj:
            if line.startswith('>'):
                name = line.split()[0][1:]
                proteins[name] = ''
            else:
                proteins[name] += line.rstrip('\n').upper()
    except:
        print('Cannot find the input file {0}. Please check the filename'.format(file))
        sys.exit(2)
    fileObj.close()
    return proteins

def writeToFile(name, file, missed, initial_peps):
    peps = initial_peps
    if missed:
        peps = missCleavage(peps, missed)
    for i in range(0, len(peps)):
        terminal = 'I'
        if i == 0:
            terminal = 'N'
        elif i == len(peps) - 1:
            terminal = 'C'
        file.write('>%s\tpeptide\t%d\tmissed=%d\t%s\t%s\n' % (key, i+1, missed, terminal, enzyme_dict[enzyme]))
        file.write('\n'.join([peps[i][j:j+60] for j in range(0,len(peps),60)]) + '\n')

if __name__ == '__main__':
    printInputInfo()
    proteins = readFile(input_file)
    writeObj = open(output_file, 'w')
    for key in proteins:
        origin = enzymeFullCut(proteins[key], enzyme_rules[enzyme])
        for i in range(0, cleavage_number+1):
            writeToFile(key, writeObj, i, origin)
    writeObj.close()