#!/usr/bin/env python3
import os
import sys
import re
import time
import pandas as pd
import numpy as np
from datetime import timedelta
from optparse import OptionParser
from collections import defaultdict
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from re import search
#Setting up script variables
parser = OptionParser()
usage = "usage: %prog [-i] csvfile [-o] outputfilename"
parser.add_option('-b', dest="binpath",help="path to the directory that only contains bins or an assembly in FASTA format", default="NA")
parser.add_option("--min_l", dest="min_l",type=float, help="minimum length of contig to write to FASTA (default = 0)", default=0)
parser.add_option("-m", dest="meta", help="If your contigs are coming from metagenomes add -m TRUE, so we can let prodigal know", default="FALSE")
parser.add_option("-o", dest="outfile", help="Basename for output files", default="NA")
parser.add_option("--bam", dest="bamfile", help="Input sorted BAM file or if from a co-assembly provide a comma-separated list of BAM files (e.g. file1.bam,file2.bam,file3.bam) (Optional layer)", default="NA")
parser.add_option("-t", dest="threads", help="Number of parallel threads used for analysis", default="1")
parser.add_option("-c", dest="CAT", help="Let us know if you want to generate taxa information for your contigs (Optional Layer) (If yes use -c TRUE)", default="FALSE")
parser.add_option("-a", dest="depthfile", help="If you already have depth file generated from jgi you can provide that file name here", default="NA")
(options, args) = parser.parse_args()
start_time_Tier_1 = time.monotonic()
#incorrect option checking
if options.binpath == "NA":
    print("You need to include the folder that contains only your paints (MAGs/assemblies) with option -b")
    raise SystemExit
if options.outfile == "NA":
    print("How will we know what to name your outputfile? Please set an output file name with option -o")
    raise SystemExit
#set in functions
#read in fasta
def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count
def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict
#calculate GC%
def GCcalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count += 1
    return (count/len(seq))*100
#calculate tetranucleotide frequency
def tetfreq(seq):
    tet_Dict = defaultdict(list)
    Nucs = ["T","A","G","C"]
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    if Nucs[a] in ["A", "G", "C", "T"] and Nucs[b] in ["A", "G", "C", "T"] and Nucs[c] in ["A", "G", "C", "T"] and Nucs[d] in ["A", "G", "C", "T"]:
                        tet = Nucs[a] + Nucs[b] + Nucs[c] + Nucs[d]
                        tet_Dict[tet] = []

    total = 0
    for e in range(len(seq)):
        f = (seq[e:e+4])
        if len(f) == 4:
            if f[0] in ["A", "G", "C", "T"] and f[1] in ["A", "G", "C", "T"] and f[2] in ["A", "G", "C", "T"] and f[3] in ["A", "G", "C", "T"]:
                tet_Dict[f].append(f)
                total += 1
    totalkmers = total
    return tet_Dict, total
#get codon usage
def codonTable(seq):
    Dict = defaultdict(lambda: defaultdict(list))
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1
    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            Dict[CodonTable[codon]][codon].append(codon)
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return Dict
#split orf names
def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length-1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]
#set dictionaries for output
Dict_summarray = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
tet_summarray = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
orf_summarray = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
taxa_summarray = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
CodonList = []
codon_summarray = defaultdict(lambda: defaultdict(lambda: 0))
#Tier 1 Upgrades
print("")
print("First Layer (sequence composition) being painted")
print("----------------------------------------------------------------------")
print("")
for bin_file in os.listdir(options.binpath):
    print("   Wokring on Bin " + bin_file)
    fasta_file = open(options.binpath + "/" + bin_file)
    fastaDict = fasta(fasta_file)
    orf = options.binpath + "/" + bin_file
    if options.meta != "FALSE":
        os.system("prodigal -a %s.%s.faa -f gff -q -o %s.%s.out -p meta -i %s -d %s.%s.nuc" % (options.outfile, bin_file, options.outfile, bin_file, orf, options.outfile, bin_file))
        os.system("rm %s.%s.out" % (options.outfile, bin_file))
    else:
        os.system("prodigal -a %s.%s.faa -f gff -q -o %s.%s.out -i %s -d %s.%s.nuc" % (options.outfile, bin_file, options.outfile, bin_file, orf, options.outfile, bin_file))
        os.system("rm %s.%s.out" % (options.outfile, bin_file))
    orf_file = open(options.outfile + "." + bin_file + ".nuc")
    orfDict = fasta(orf_file)
    orfCountDict = defaultdict(list)
    for key in fastaDict.keys():
        seq = fastaDict[key]
        seq_len = len(seq)
        if seq_len <= options.min_l:
            pass
        else:
            gc_per = GCcalc(seq)
            Dict_summarray[key]["gc_per"] = gc_per
            Dict_summarray[key]["bin"] = bin_file
            Dict_summarray[key]["contig_length"] = seq_len
            #Dict_summarray[key]["orf_count"] = len(orfCountDict[key])
            tetfreqAndtotalkmers = tetfreq(seq)
            tet_freq = tetfreqAndtotalkmers[0]
            for h in tet_freq.keys():
                tet_summarray[key][h] = len(tet_freq[h])/tetfreqAndtotalkmers[1]

    for f in orfDict.keys():
        contig = allButTheLast(f,"_")
        if contig in Dict_summarray.keys():
            orfCountDict[contig].append(f)

    for a in orfDict.keys():
        contig = allButTheLast(a, "_")
        if contig in Dict_summarray.keys():
            seq = (orfDict[a])
            ct = (codonTable(seq))
            for b in ct.keys():
                for c in ct[b]:
                    orf_summarray[contig][b][c].append(len(ct[b][c]))

    for d in sorted(orf_summarray.keys()):
        contig = d
        if contig in Dict_summarray.keys():
            for e in sorted(orf_summarray[d]):
                localtotal = 0
                for f in sorted(orf_summarray[d][e]):
                    co = sum(orf_summarray[d][e][f])
                    localtotal += co

                for g in sorted(orf_summarray[d][e]):
                    co = sum(orf_summarray[d][e][g])
                    codon_summarray[d][g] = co / localtotal
                    if g not in CodonList:
                        CodonList.append(g)

print("")
print("Computing Principal component analysis on Tetranucleotide Frequency")
print("note : This is 4kmer_PCA1, 4kmer_PCA2, and 4kmer_PCA3 in the tsv file ")
print("")
#4kmer PCA
tet_df = pd.DataFrame.from_dict(tet_summarray,orient='index')
tet_df_w_labels=tet_df.rename_axis("contig").reset_index()
tet_df_standard = StandardScaler().fit_transform(tet_df)
tet_pca = PCA(n_components=3)
tet_principalComponents = tet_pca.fit_transform(tet_df_standard)
tet_principalDf = pd.DataFrame(data = tet_principalComponents,columns = ['4kmer_PCA1', '4kmer_PCA2', '4kmer_PCA3'])
tet_pca = pd.concat([tet_principalDf, tet_df_w_labels], axis=1)
tet_pca = tet_pca[tet_pca.columns[0:4]]
tet_pca = tet_pca.set_index('contig')
tet_pca.index.names = [None]
tet_pca_summarray = tet_pca.to_dict(orient="index")
print("")
print("Computing Principal component analysis on Codon Usage Bias")
print("note : This is codon_PCA1, codon_PCA2, and codon_PCA3 in the tsv file")
print("")
#codon bias PCA
codon_df = pd.DataFrame.from_dict(codon_summarray,orient='index')
codon_df_no_NaN = codon_df.fillna(0)
codon_df_w_labels=codon_df_no_NaN.rename_axis("contig").reset_index()
codon_df_standard = StandardScaler().fit_transform(codon_df_no_NaN)
codon_pca = PCA(n_components=3)
codon_principalComponents = codon_pca.fit_transform(codon_df_standard)
codon_principalDf = pd.DataFrame(data = codon_principalComponents,columns = ['codon_PCA1', 'codon_PCA2', 'codon_PCA3'])
codon_pca = pd.concat([codon_principalDf, codon_df_w_labels], axis=1)
codon_pca = codon_pca[codon_pca.columns[0:4]]
codon_pca = codon_pca.set_index('contig')
codon_pca.index.names = [None]
codon_pca_summarray = codon_pca.to_dict(orient="index")
print("")
print("Computing Principal component analysis on Tier 1 Sequence Composition Metrics")
print("note : This is x,y,z in the tsv file (x and y are automatically plotted)")
print("")
#all metric PCA
GC_df = pd.DataFrame.from_dict(Dict_summarray,orient='index')
GC_df_w_labels=GC_df.rename_axis("contig").reset_index()
pca_ALL = pd.concat([codon_df_no_NaN,tet_df,GC_df], axis=1)
pca_ALL = pca_ALL[pca_ALL.columns[0:321]]
pca_ALL_w_labels=pca_ALL.rename_axis("contig").reset_index()
pca_ALL_standard = StandardScaler().fit_transform(pca_ALL)
ALL_pca = PCA(n_components=3)
ALL_principalComponents = ALL_pca.fit_transform(pca_ALL_standard)
ALL_principalDf = pd.DataFrame(data = ALL_principalComponents,columns = ['x', 'y', 'z'])
ALL_pca = pd.concat([ALL_principalDf, pca_ALL_w_labels], axis=1)
ALL_pca = ALL_pca[ALL_pca.columns[0:4]]
ALL_pca = ALL_pca.set_index('contig')
ALL_pca.index.names = [None]
ALL_pca_summarray = codon_pca.to_dict(orient="index")
end_time_Tier1 = time.monotonic()
print("Finished in " + str(timedelta(seconds=end_time_Tier1 - start_time_Tier_1)))
#calculating depth
if options.bamfile != "NA":
    start_time_Tier_2 = time.monotonic()
    print("")
    print("Second Layer (depth) being painted")
    print("----------------------------------------------------------")
    print("")
    bamfiles = options.bamfile
    bamfileList = bamfiles.split(",")
    bamString = ''
    #print(bamfileList2)
    for bam in bamfileList:
        bamString += bam
        bamString += " "
    os.system("jgi_summarize_bam_contig_depths --outputDepth %s.depth %s --minContigLength %s" % (options.outfile, bamString, options.min_l))
    end_time_Tier2 = time.monotonic()
    print("Finished in " + str(timedelta(seconds=end_time_Tier2 - start_time_Tier_2)))

#get taxonomic information on contigs
if options.CAT != "FALSE":
    start_time_Tier_3 = time.monotonic()
    print("")
    print("Third Layer (taxonomy) being painted")
    print("----------------------------------------------------------")
    print("")
    CATdb = os.path.abspath(os.path.dirname(sys.argv[0])) + "/CAT_prepare_20210107"
    rscript = os.path.abspath(os.path.dirname(sys.argv[0])) + "/third_layer.R"
    os.system("cat %s/* > %s.all.contigs.fa" % (options.binpath, options.outfile))
    print(" Classifiying contigs using CAT")
    print(" Please make sure the you cite them when using this program:")
    print(" von Meijenfeldt FAB, Arkhipova K, Cambuy DD, Coutinho FH, Dutilh BE.")
    print(" Robust taxonomic classification of uncharted microbial sequences and bins with CAT and BAT.")
    print(" Genome Biology. 2019;20:217.  https://doi.org/10.1186/s13059-019-1817-x")
    print("")
    os.system("CAT contigs -c %s.all.contigs.fa -d %s/2021-01-07_CAT_database/ -t %s/2021-01-07_taxonomy/ -o %s.CAT_out -n %s --no_stars --quiet --top 10 --I_know_what_Im_doing" % (options.outfile, CATdb, CATdb, options.outfile, options.threads))
    os.system("CAT add_names -i %s.CAT_out.contig2classification.txt -o %s.classified_contigs -t %s/2021-01-07_taxonomy/ --exclude_scores -q --only_official --quiet" % (options.outfile, options.outfile, CATdb))
    os.system("Rscript %s %s.classified_contigs" % (rscript, options.outfile))
    #######################################
    taxa = open("%s.classified_contigs" % options.outfile)
    for g in taxa:
       h = g.rstrip().split("\t")
       taxa_summarray[h[0]]["Taxa_Rank_1"] = h[1]
       taxa_summarray[h[0]]["Taxa_Rank_2"] = h[2]
       taxa_summarray[h[0]]["Taxa_Rank_3"] = h[3]
       taxa_summarray[h[0]]["Taxa_Rank_4"] = h[4]
       taxa_summarray[h[0]]["Taxa_Rank_5"] = h[5]
       taxa_summarray[h[0]]["Taxa_Rank_6"] = h[6]
       taxa_summarray[h[0]]["Taxa_Rank_7"] = h[7]
       taxa_summarray[h[0]]["Taxa_Rank_8"] = h[8]
       taxa_summarray[h[0]]["Taxa_Rank_9"] = h[9]
       taxa_summarray[h[0]]["Taxa_Rank_10"] = h[10]
       taxa_summarray[h[0]]["Taxa_Rank_11"] = h[11]
    end_time_Tier3 = time.monotonic()
    print("Finished in " + str(timedelta(seconds=end_time_Tier3 - start_time_Tier_3)))

#Writing output files for barbizon input
print("")
print("Printing output files for your enjoyment!")
print("------------------------------------------------------")
print("")
out = open(options.outfile + ".barbizon_input.tsv", "w")
out.write("contigName\tBin\tLength\tGC\t4kmer_PCA1\t4kmer_PCA2\t4kmer_PCA3\tcodon_PCA1\tcodon_PCA2\tcodon_PCA3\tx\ty\tz\n")
for j in Dict_summarray.keys():
    out.write(j + "\t" + Dict_summarray[j]["bin"] + "\t" + str(Dict_summarray[j]["contig_length"]) + "\t" + str(round(Dict_summarray[j]["gc_per"], 4)))
    for k in tet_pca_summarray[j]:
        out.write("\t" + str(tet_pca_summarray[j][k]))
    for l in codon_pca_summarray[j]:
        out.write("\t" + str(codon_pca_summarray[j][l]))
    for m in ALL_pca_summarray[j]:
        out.write("\t" + str(codon_pca_summarray[j][m]))
    out.write("\n")
out.close()
os.system("rm %s.*.faa %s.*.nuc" % (options.outfile, options.outfile))

if options.bamfile != "NA":
    barbizon_input = options.outfile + ".barbizon_input.tsv"
    temp_input = pd.read_csv(barbizon_input, sep = "\t")
    depth = options.outfile + ".depth"
    temp_depth = pd.read_csv(depth, sep = "\t")
    output_depth = pd.merge(temp_input,temp_depth, on='contigName')
    output_depth = output_depth[output_depth.columns.drop(list(output_depth.filter(regex='-var')))]
    output_depth = output_depth[output_depth.columns.drop(list(output_depth.filter(regex='contigLen')))]
    output_depth = output_depth[output_depth.columns.drop(list(output_depth.filter(regex='totalAvgDepth')))]
    output_depth = output_depth.dropna()
    depth_output_name = options.outfile + ".barbizon_input_wCov.tsv"
    output_depth.to_csv(depth_output_name, sep="\t")
    os.system("mv %s.barbizon_input_wCov.tsv %s.barbizon_input.tsv" % (options.outfile, options.outfile))
    os.system("rm %s.depth" % (options.outfile))

if options.depthfile != "NA":
    barbizon_input = options.outfile + ".barbizon_input.tsv"
    temp_input = pd.read_csv(barbizon_input, sep="\t")
    depth = pd.read_csv(options.depthfile, sep="\t")
    output_depth = pd.merge(temp_input, depth, on='contigName')
    output_depth = output_depth[output_depth.columns.drop(list(output_depth.filter(regex='-var')))]
    output_depth = output_depth[output_depth.columns.drop(list(output_depth.filter(regex='contigLen')))]
    output_depth = output_depth[output_depth.columns.drop(list(output_depth.filter(regex='totalAvgDepth')))]
    output_depth = output_depth.dropna()
    depth_output_name = options.outfile + ".barbizon_input_wCov.tsv"
    output_depth.to_csv(depth_output_name, sep="\t")
    os.system("mv %s.barbizon_input_wCov.tsv %s.barbizon_input.tsv" % (options.outfile, options.outfile))

if options.CAT != "FALSE":
    inputTSV = open(options.outfile + ".barbizon_input.tsv", "r")
    out = open(options.outfile + ".barbizon_input_wCAT.tsv", "w")
    for i in inputTSV:
        ls = i.rstrip().split("\t")
        if ls[0] == "id":
            out.write(i.rstrip() + "\ttaxon|d\tTaxa_Rank_2\tTaxa_Rank_3\tTaxa_Rank_4\tTaxa_Rank_5\tTaxa_Rank_6\tTaxa_Rank_7\tTaxa_Rank_8\tTaxa_Rank_9\tTaxa_Rank_10\tTaxa_Rank_11\n")
        else:
            out.write(i.rstrip() + "\t" + taxa_summarray[ls[0]]["Taxa_Rank_1"] + "\t" + taxa_summarray[ls[0]]["Taxa_Rank_2"] + "\t"+ taxa_summarray[ls[0]]["Taxa_Rank_3"] + "\t" +
                      taxa_summarray[ls[0]]["Taxa_Rank_4"] + "\t"+ taxa_summarray[ls[0]]["Taxa_Rank_5"] + "\t"+ taxa_summarray[ls[0]]["Taxa_Rank_6"] + "\t"+ taxa_summarray[ls[0]]["Taxa_Rank_7"] + "\t" +
                      taxa_summarray[ls[0]]["Taxa_Rank_8"] + "\t"+ taxa_summarray[ls[0]]["Taxa_Rank_9"] + "\t" + taxa_summarray[ls[0]]["Taxa_Rank_10"] + taxa_summarray[ls[0]]["Taxa_Rank_11"] + "\n")
    out.close()
    os.system("mv %s.barbizon_input_wCAT.tsv %s.barbizon_input.tsv" % (options.outfile, options.outfile))
    #os.system("rm %s.CAT*" % (options.outfile))
