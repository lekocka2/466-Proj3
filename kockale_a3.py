#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 09:47:49 2021

@author: logankocka
"""
# project 3
import sys
import csv
import pandas
import seaborn
from pylab import savefig

def get_seq_fasta(file):
    # name_list=[]; desc_list=[]; seq_list=[]
    dict1 = {}
    # print("A Inside readFile() [{0}]".format(file))
    fh = open(file, 'r')
    name, desc, seq = ('', '', '')
    for line in fh:
        line = line.strip()

        if line.startswith('>'):
            if seq !='':
                name, desc, seq = ('', '', '')
            name = line[1:line.index(' ')]
            desc = line[line.index(' ')+1:]
        else:
            seq=seq+line
            dict1[name] = seq # key is the name, value is the sequence
            
    return (dict1)


def get_expression(fileName1, fileName2):
    dict_lung, dict_prot = ({}, {}) # keys are gene names, values are differentials
    
    for gene in csv.reader(open(fileName1), delimiter='\t'):
        if gene[0].startswith("E"):
            if int(gene[1])/int(gene[2]) >= 1.5: 
                temp = "+"
            elif int(gene[1])/int(gene[2]) <= 0.667: 
                temp = "-"
            else: 
                temp = "."
        
            l = gene[1], ":", gene[2], ":", temp
            dict_lung[gene[0]] = "".join(l)

    for gene in csv.reader(open(fileName2), delimiter='\t'):
        if gene[0].startswith("E"):
            if int(gene[1]) > int(gene[2]): 
                temp = "+"
            elif int(gene[1]) == int(gene[2]): 
                temp = "."
            else: 
                temp = "-"
        
            l = gene[1], ":", gene[2], ":", temp
            dict_prot[gene[0]] = "".join(l)
    
    return dict_lung, dict_prot


def common_genes(dict1, dict2):
    set1 = set(dict1.keys())
    set2 = set(dict2.keys())
    set3 = set1.intersection(set2)
    common_set = sorted(set3)
    
    return common_set


def unique_genes(dict1, dict2):
    common_set = common_genes(dict1, dict2)
    unique_lung = set(dict1.keys())
    unique_prot = set(dict2.keys())
    # remove items in common set
    for i in common_set:
        unique_lung.remove(i)
        unique_prot.remove(i)
    
    return sorted(unique_lung), sorted(unique_prot)


def compare(dict1, dict2, common_set):
    for item in common_set:
        if dict1[item][-1] == dict2[item][-1]:
            temp = "Congruence"
        else: temp = "Disparity"
        
        print(item, " NumL:NumC:Exp (", dict1[item], ") NumP:NumC:Exp (", dict2[item], ") ", temp, sep="")

    return


def get_nuc_pos_profile(dictionar, gene_set, start, end): 
    # print positions on one line
    # loop through the genes in the set
    # tally nucleotides
    # listPos = []
    listPos = range(start, end+1)
    # listPos = [str(i) for i in listPos]
    print("Position ", end="")
    # for num in listPos:
    print("".join([str(i).center(4) for i in listPos]))
    print("--------  ", end = "")
    # for list of dashes
    length = end-start+1
    listDash = []
    for i in range(0, length):
        listDash.append('---')
    # print(listDash)
    print(" ".join([str(i).center(2) for i in listDash]))
    
    # # use lists to tally
    A = []; T = []; G = []; C = []; # initialize
    A_count = G_count = C_count = T_count = 0 # initialize
   
    # try with dataframes
    temp_list = []
    new_dict = {}
    for gene in gene_set: # this loops thru the genes in the gene list
        for char in dictionar[gene]: # this loops thru the chars in a sequence
            temp_list.append(char)
        new_dict[gene] = temp_list
        temp_list = []
    
    values_only = new_dict.values() # get list of values for dataframe
    df = pandas.DataFrame(values_only)
    # now I have a data frame with the characters of each sequence in rows with position as column headers
    # print("DF", df)
    # now count down columns each char counts
    for col in df.columns[start-1:end]:
        # A_count = df[col].value_counts()["A"]
        A_count = len(df[df[col] == 'A'])
        A.append(A_count)
        T_count = len(df[df[col] == 'T'])
        T.append(T_count)
        G_count = len(df[df[col] == 'G'])
        G.append(G_count)
        C_count = len(df[df[col] == 'C'])
        C.append(C_count)
        
        A_count = G_count = C_count = T_count = 0  # reset the values to zero
        
    # print("A", A)
    # print("T", T) # seems to work 
    print("A         ", end="")
    print("  ".join([str(i).center(2) for i in A]))
    print("T         ", end="")
    print("  ".join([str(i).center(2) for i in T]))
    print("G         ", end="")
    print("  ".join([str(i).center(2) for i in G]))
    print("C         ", end="")
    print("  ".join([str(i).center(2) for i in C]))
    print()
    return


def heatmap (df1, df2, imgfile):  
    # hm = heatmap()
    df1["fold_L"] = df1.NumL / df1.NumC
    df2["fold_P"] = df2.NumP / df2.NumC
    # this is the one I plot from
    df_plot = pandas.merge(df1, df2, how='inner', on=['GeneID'])
    # print(df_plot)
    rowNames = df_plot["GeneID"] # save this column for row names
    df_plot = df_plot.drop(df_plot.columns[[0,1,2,4,5]], axis=1)
    # insert row names
    df_plot.index = rowNames
    
    fig = seaborn.heatmap(df_plot, annot=True)
    fig.set_title("Lung and Prostate Cancer Expression")
    fig.figure.subplots_adjust(left = 0.3)
    figure = fig.get_figure()    
    figure.savefig(imgfile, dpi=400)
    return


def main():
    args = sys.argv
    # user = args[0].split("_")[0]
    fileName = args[1]
    txtFile1 = args[2]
    txtFile2 = args[3]
    
    print("1. Files: FASTA ({0}) Cancer Expression ({1}, {2})".format(fileName, txtFile1, txtFile2))
    dictionar = get_seq_fasta(fileName)
    
    print("\n2. Differentially expressed gene detected for lung cancer-control comparison:")
    lung, prost = get_expression(txtFile1, txtFile2)
    print("GeneID     NumL:NumC:Exp")
    for key, value in lung.items():
        print(key, value)

    print("\n3. Differentially expressed gene detected for prostate cancer-control comparison:")
    print("GeneID     NumP:NumC:Exp")
    for key, value in prost.items():
        print(key, value)

    print("\n4. Expression comparison (lung vs control) vs (prostate vs control)")
    compare(lung, prost, common_genes(lung, prost))

    result = common_genes(lung, prost)
    print("\n5. Common DEGs (expressed genes) for both lung and prostate cancer\n", ", ".join(list(result)))

    uniqueL, uniqueP = unique_genes(lung, prost)
    print("\n6. DEGs uniquely expressed in lung cancer \n", ", ".join(list(uniqueL)))

    print("\n7. DEGs uniquely expressed in prostate cancer \n", ", ".join(list(uniqueP)))
    
    examChoice = input("\n8. Which gene set do you want to examine? (5, 6, or 7) ")
    print("\n8. Your answer is [{0}]".format(examChoice))
    if examChoice == "5":
        print("\nWe will extract nucleotide position profile for the following sequences:")
        print(", ".join(list(result)))
        use = result
    elif examChoice == "6":
        print("\nWe will extract nucleotide position profile for the following sequences:")
        print(", ".join(list(uniqueL)))
        use = uniqueL
    elif examChoice == "7":
        print("\nWe will extract nucleotide position profile for the following sequences:")
        print(", ".join(list(uniqueP)))
        use = uniqueP
    elif examChoice == "":
        print("\nWe will extract nucleotide position profile for the following sequences:")
        set1 = set(lung.keys())
        set2 = set(prost.keys())
        All = set1.union(set2)
        print(", ".join(list(All)))
        use = All
    else:
        print("Invalid choice.")
        sys.exit()
    len = input("\nPlease enter the start and end potions (e.g., 1-140) ")
    start = int(len.split("-")[0])
    end = int(len.split("-")[1])
    if start < 1 | end > 140:
        print("Invalid entry.")
        sys.exit()
    else:
        get_nuc_pos_profile(dictionar, use, start, end)
    
    
    # heatmap codes -------------------------
    df1 = pandas.read_csv(txtFile1, sep="\t")
    df2 = pandas.read_csv(txtFile2, sep="\t")
    imgfile = 'Lung_prostate_heatmap'  # imgfile = ???
    heatmap(df1, df2, imgfile)
    
    return


if __name__=='__main__':
    main()
    