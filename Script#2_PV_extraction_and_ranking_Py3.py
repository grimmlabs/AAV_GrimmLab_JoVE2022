#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 18:47:02 2022

@author: kleopatrarapti
"""

import os
import pandas as pd
from Bio.Seq import Seq
from datetime import datetime

config={}
exec(open("//mypath/NGS_Script/Barcode_Script_JoVE.conf").read())
print(my_dir)
print("-------------------------------------Script#2 is running-------------------------------------")

with open (filename_sample_file) as temp:
   sample=dict(line.strip().split() for line in temp if line.strip())


objects = [file for file in os.listdir(my_dir) if file.endswith('.txt')]
for filename in objects:
    print("file being analyzed:", 2*"\t", filename, 2*"\t", datetime.now(), "\n")
#at df one can also choose the number of rows to analyze nrows=5000
    df = pd.read_csv(my_dir+filename, skiprows=5, header = 0, sep = '\t')
    dfhead = pd.read_csv(my_dir+filename, skiprows=1, header=None,index_col=None, nrows=2, sep = '\t')
    PV_DNAs = list(df['Sample:'])
#make a dictionary to get pair of DNA and read number values. This will be used to count the valid reads
    dict_DNA_Reads = dict(zip(df['Sample:'], df['# ']))

    PV_PRs =[]
    PV_frwORrev = []
    RG = "AGAGGC"
    SG = "AGTGGC"

    valid_reads=0
    invalid_reads=0
    for PV_DNA in PV_DNAs:
        PV_DNA_rev = str(Seq(PV_DNA).reverse_complement())
#        if (PV_DNA[:6].find(RG)!=-1) or (PV_DNA[:6].find(SG)!=-1):
        if (PV_DNA[:6]==RG) or (PV_DNA[:6]==SG):
#            df.loc['# ', df[df['Pid'] == 'p01']]
            valid_reads+=dict_DNA_Reads[PV_DNA]
            PV_PRs.append(str(Seq(PV_DNA).translate()))
            PV_frwORrev.append("frw")
        elif (PV_DNA_rev[:6]==RG) or (PV_DNA_rev[:6]==SG):
            valid_reads+=dict_DNA_Reads[PV_DNA]
            PV_PRs.append(str(Seq(PV_DNA_rev).translate()))
            PV_frwORrev.append("rev")
        else:
            PV_PRs.append("not valid")
            PV_frwORrev.append("NA")
            invalid_reads+=1
    df["Frw or Rev"]=PV_frwORrev
    df["PVs"]=PV_PRs

#df dataframe now contains 4 columns (Sample:,# ,Frw or Rev,PVs) and all the original DNA sequences. The  
#Generate the dfvalid dataframe with just valid reads by selected all but "not valid" from column PVs/4 columns
    dfvalid = df[df["PVs"]!="not valid"]
#Generate the dfPVs dataframe from dfvalid with unique PVs/2 columns
    dfPVs = dfvalid[['# ', 'PVs']].groupby('PVs')['# '].agg(['sum','count']).sort_values(by=['sum'],ascending=False).reset_index(level=0)
    dfPVs.rename(columns = {'sum':'# '}, inplace = True)

#Generate a table with the statistics called dfhead. Add to the statistics from the original file which are in dfhead
    row3 = pd.Series(["# of Valid PV reads: ", str(valid_reads)])
    row4 = pd.Series(["# of Invalid PV reads: ", str(invalid_reads)])
    #remove 1 count from the unique values of PV_PRs to count for "not valid"
    row5 = pd.Series(["# of unique PV reads: ", str(len(set(PV_PRs))-1)])                  
    dfheadnew=dfhead.append([row3,row4,row5],ignore_index=True)
    
#Generate a file with just the PVs that are unique
    dfheadnew[2]=float("NaN")
    col_names=df.columns[:2]
    dfheadnew.columns=[df.columns[3], df.columns[1],'count']
    dfPVsexport = dfheadnew.append(dfPVs.sort_values(by=['# '],ascending=False), sort=False, ignore_index =False)

#Generate the rest of the files
    dfheadnew[3]=float("NaN")
    dfheadnew.columns=df.columns
#Generate the file with all data
    dfallexport = dfheadnew.append(df.sort_values(by=['# '],ascending=False), sort=False, ignore_index =False)
#Generate the file with valid data
    dfvalidexport = dfheadnew.append(dfvalid.sort_values(by=['# '],ascending=False), sort=False, ignore_index =False)

#replace filename with new name from Zuordnung.txt
    if filename in sample.keys():
        newname = sample[filename]

    dfPVsexport.to_excel(my_dir + newname + "_PVs.xlsx")
    dfallexport.to_excel(str(my_dir)+newname+"_all.xlsx")
    dfvalidexport.to_excel(str(my_dir)+newname+"_validSeq.xlsx")

