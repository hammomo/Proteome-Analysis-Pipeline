#!/usr/bin/python
# -*- coding: UTF-8 -*-

m_map = {'A' :  71.0371, 'C' : 103.0092, 'D' : 115.0269, 'E' : 129.0426,
        'F' : 147.0684, 'G' :  57.0215, 'H' : 137.0589, 'I' : 113.0841,
        'K' : 128.0950, 'L' : 113.0841, 'M' : 131.0405, 'N' : 114.0429,
        'P' :  97.0528, 'Q' : 128.0586, 'R' : 156.1011, 'S' :  87.0320,
        'T' : 101.0477, 'V' :  99.0684, 'W' : 186.0793, 'Y' : 163.0633,
        '\s' : 0.0, '*' : 0.0
       }
a_map = {'A' :  71.08, 'C' : 103.14, 'D' : 115.09, 'E' : 129.12,
        'F' : 147.18, 'G' :  57.05, 'H' : 137.14, 'I' : 113.16,
        'K' : 128.17, 'L' : 113.16, 'M' : 131.19, 'N' : 114.10,
        'P' :  97.12, 'Q' : 128.13, 'R' : 156.19, 'S' :  87.08,
        'T' : 101.10, 'V' :  99.13, 'W' : 186.21, 'Y' : 163.18,
        '\s' : 0.0, '*' : 0.0
       }

water_m_mass = 18.0106
water_a_mass = 18.0153

def monoisotopicMass(seq):
       sum = 0
       for aa in seq:
              sum +=  m_map[aa]
       sum += water_m_mass
       return sum

def averageMass(seq):
       sum = 0
       for aa in seq:
              sum += a_map[aa]
       sum += water_a_mass
       return sum

if __name__ == '__main__':
       readObj = open('peptides.peps', 'r')
       monoObj = open('monoisotopic-mass.pepmasses', 'w')
       avgObj = open('average-mass.pepmasses', 'w')
       p = readObj.readline().split(':')[1].rstrip('\n')
       line = readObj.readline()
       port_name = ''
       peptide = ''
       sequence = ''
       z = 1
       while line:
              if line.startswith('>'):
                     port_name = line.split()[0].split('_')[0].replace('>', '')
                     peptide = line.split()[-1].rstrip('\n')
              else:
                     sequence = line.rstrip('\n')
                     m_mass = monoisotopicMass(sequence)
                     a_mass = averageMass(sequence)
                     monoObj.write('{0}\t{1:>10s}\t{2:>10.4f}\t{3:>2d}\t{4:>1s}\t{5}\n'.format(port_name, peptide, m_mass, z, p, sequence))
                     avgObj.write('{0}\t{1:>10s}\t{2:>10.4f}\t{3:>2d}\t{4:>1s}\t{5}\n'.format(port_name, peptide, a_mass, z, p, sequence))
              line = readObj.readline()
       readObj.close()
       monoObj.close()
       avgObj.close()