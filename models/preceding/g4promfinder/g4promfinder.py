#Source: https://github.com/MarcoDiSalvo90/G4PromFinder/blob/master/G4PromFinder.py

import numpy as np
import matplotlib.pyplot as plt
import re
from Bio import SeqIO
from Bio.Seq import Seq
from openpyxl import load_workbook
from Bio import motifs
MC = []

print("G4PromFinder is a genome-wide predictor of transcription promoters in bacterial genomes. It analyzes both the strands. It is recommended for GC-rich genomes\n\n")
print("Input: Genome file")
print("Output: a text file containing promoter coordinates and a file containing informations about them")
print("Genome file must be in fasta format")
File_genome = input("Enter the input genome file name: ")

try:
    for seq_record in SeqIO.parse(File_genome, "fasta"):
        MC.append(str(seq_record.seq))
    print ("G4PromFinder is working, please wait...")
    genome = MC[0] 
    genome = genome.upper()
    gc_total = 100*(genome.count("G")+genome.count("C"))/len(genome)
    pattern1 = "G{2,4}\D{1,10}G{2,4}\D{1,10}G{2,4}\D{1,10}G{2,4}"
    pattern2 = "C{2,4}\D{1,10}C{2,4}\D{1,10}C{2,4}\D{1,10}C{2,4}"
    patternA = "TA\D{3}T"
    patternB = "A\D{3}TA"
    patternC = "TTGAC"
    patternD = "GTCAA"
    predictions = 0
    cod1 = 1
    cod2 = 1
    AT_content = []
    GC_content = []
    SI1 = 0
    SI2 = 0
    NO1 = 0
    NO2 = 0
    q = open("Promoter coordinates.txt","a")
    q.write("Region           Strand              Start               End\n")
    
    scale = 50
    for j in range(100000000):
        x1 = scale + j
        x2 = scale + j + 25
        OK = 0
        if x2 > len(genome):
            break
        else:
            w = genome[x1:x2]
            at =100*(w.count("A")+w.count("T"))/len(w)
            prom = genome[x1-50:x2]
            if at >= 40:
                if re.findall(pattern1,prom) or re.findall(pattern2,prom):
                    if re.findall(pattern1,prom):
                        ricerca = re.findall(pattern1,prom)
                        for lk in ricerca:
                            if len(lk) <= 30:
                                OK += 1
                                break
                    elif re.findall(pattern2,prom):
                        ricerca = re.findall(pattern2,prom)
                        for lk in ricerca:
                            if len(lk) <= 30:
                                OK += 1
                                break
                    if OK == 1:
                        predictions += 1
                        dati = []
                        for g in range(50):
                            if x2+g <= len(genome):
                                w = genome[x1+g:x2+g]
                                at =100*(w.count("A")+w.count("T"))/len(w)
                                prom = genome[x1+g-50:x2+g]
                                if at >= 40:
                                    if re.search(pattern2,prom) or re.search(pattern1,prom):
                                        dati.append(at)
                                    else:
                                        dati.append(0)
                                else:
                                    dati.append(0)
                            else:
                                break
                        maxP = np.argmax(dati)
                        at = np.max(dati)
                        x1 = x1 + maxP
                        x2 = x2 + maxP
                        scale = x2 - j - 1 + 50
                        prom = genome[x1-50:x2]
                        gc = 100*(prom.count("G")+prom.count("C"))/len(prom)
                        q.write("\nP")
                        q.write(str(cod1))
                        q.write("plus             ")
                        q.write("positivo           ")
                        q.write(str(x1-50))
                        q.write("           ")
                        q.write(str(x2))
                        cod1 += 1
                        AT_content.append(at)
                        GC_content.append(gc)
                        if re.search(patternA,prom):
                            SI1 += 1
                        else:
                            NO1 += 1
                        if re.search(patternC,prom):
                            SI2 += 1
                        else:
                            NO2 += 1
                            
    scale = 0
    for j in range(100000000):
        x1 = scale + j
        x2 = scale + j + 25
        OK = 0
        if x2 > len(genome)-50:
            break
        else:
            w = genome[x1:x2]
            at =100*(w.count("A")+w.count("T"))/len(w)
            prom = genome[x1:x2+50]
            if at >= 40:
                if re.findall(pattern1,prom) or re.findall(pattern2,prom):
                    if re.findall(pattern1,prom):
                        ricerca = re.findall(pattern1,prom)
                        for lk in ricerca:
                            if len(lk) <= 30:
                                OK += 1
                                break
                    elif re.findall(pattern2,prom):
                        ricerca = re.findall(pattern2,prom)
                        for lk in ricerca:
                            if len(lk) <= 30:
                                OK += 1
                                break
                    if OK == 1:
                        dati = []
                        predictions += 1
                        for g in range(50):
                            if x2+g <= len(genome):
                                w = genome[x1+g:x2+g]
                                at =100*(w.count("A")+w.count("T"))/len(w)
                                prom = genome[x1+g:x2+g+50]
                                if at >= 40:
                                    if re.search(pattern1,prom) or re.search(pattern2,prom):
                                        dati.append(at)
                                    else:
                                        dati.append(0)
                                else:
                                    dati.append(0)
                            else:
                                break
                        maxP = np.argmax(dati)
                        at = np.max(dati)
                        x1 = x1 + maxP
                        x2 = x2 + maxP 
                        scale = x2 - j - 1 + 50
                        prom = genome[x1:x2+50]
                        gc = 100*(prom.count("G")+prom.count("C"))/len(prom)
                        q.write("\nP")
                        q.write(str(cod2))
                        q.write("minus             ")
                        q.write("negativo           ")
                        q.write(str(x1))
                        q.write("           ")
                        q.write(str(x2+50))
                        cod2 += 1
                        AT_content.append(at)
                        GC_content.append(gc)
                        if re.search(patternB,prom):
                            SI1 += 1
                        else:
                            NO1 += 1
                        if re.search(patternD,prom):
                            SI2 += 1
                        else:
                            NO2 += 1
                            
    q.close()
    
    with open("ABOUT PROMOTERS.txt", "a") as file:
        file.write("\n\nMean GC% content of genome sequence: ")
        file.write(str(gc_total))
        file.write("%")
        file.write("\n\nTotal number of prediction in the analyzed sequence: ")
        file.write(str(predictions))
        file.write("\n\n\nMean GC-content of putative promoters: ")
        file.write(str(np.mean(GC_content)))
        file.write("%\nMean AT-content of the putative promoters AT-rich elements: ")
        file.write(str(np.mean(AT_content)))
        file.write("%\n\n\nPercent of putative promoters with the -10 consensus TANNNT: ")
        file.write(str(100*SI1/(SI1+NO1)))
        file.write("%\nPercent of putative promoters with the -35 consensus TTGAC: ")
        file.write(str(100*SI2/(SI2+NO2)))
        file.write("%\n\n\n")
    
    print("Work finished, see output files in the current directory")
    
except IOError:   
    print ("File %s inexistent in the current directory!" %(File_genome))
                        