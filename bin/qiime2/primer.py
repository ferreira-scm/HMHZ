#!/usr/bin/env python3

from Bio.Seq import Seq
import os

path="/SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/amplicon/"
path2="/SAN/Susanas_den/gitProj/HMHZ/tmp/qiime2/primer/"
#os.mkdir(path)
#os.mkdir(path2)

with open("/SAN/Susanas_den/gitProj/HMHZ/data/primer_list.csv", "r") as primerlist:
    next(primerlist)
    for line in primerlist:
        line = line.rstrip()
        primerF,primerR,seqF,seqR=line.split(",")
        os.mkdir(path + primerF + primerR)
        seqF=Seq(seqF)
        seqR=Seq(seqR)
        RRC=seqR.reverse_complement()
        FRC=seqF.reverse_complement()
        a=("-a %s...%s" % (seqF,RRC))
        A=("-A %s...%s" % (seqR, FRC))
        f= open(path2+primerF+primerR, "w")
        f.write(a + " " + A)
        f.close
        
#save each line is a text file in a  new folder
