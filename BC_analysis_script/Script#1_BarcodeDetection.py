#!/usr/bin/env python
#__________________________________DNA-BARCODE-SCRIPT_____________________________________________ 
#
#The script performes identification of the barcode sequences by flanking constant hexamers (sequence for hexamers can be modified in configuration file)  
import gzip
import os
import sys

print("\n"*5+"====== Barcode Script 1.4 (modified by J. Sippel, J. Weinmann, S. Weis and O. Maiakovska) ======\n\n")

#open config-file with variables 
config={}
exec(open(sys.argv[1]) .read())
samples={"":""}
print(vars())
j=0
#if samples_name exists, the samples and their corresponding index sequences are stored in a dictionary called samples: {index:sample}; if sample_name is empty, all values are set to 0/''
#samples_name -> demult_indexes
if not "samples_name" in vars() or samples_name=="":
    samples_name=""
    BCS1_size=0
    BCS2_size=0
    BCS1_right=""
    BCS2_left=""  
else:  
    with open(samples_index_file) as temp:
        samples=dict(line.strip().split() for line in temp if line.strip())

#variants and their corresponding barcode-sequences are stored in a dictionary called variants: {'barcode':'variant'}
with open(variants_barcode_file) as temp:
    print("variants_barcode_file = " + variants_barcode_file)
    variants=dict(line.strip().split() for line in temp if line.strip())
    print(len(list(variants.keys())[0]))
    
#contaminations and their barcode-sequences are stored in a dictionary called contaminationsraw: {'barcode':'contamination'}
with open(contaminations_barcode_file) as temp:
    print("contaminations_barcode_file = " + contaminations_barcode_file)
    contaminationsraw = dict(line.strip().split() for line in temp if line.strip())
    
#remove entries in {contaminationsraw} that are identical with applied barcodes in variants
contaminations = {}
for c in contaminationsraw:
    if c not in variants:
        d = {c:contaminationsraw[c]}
        contaminations.update(d)

#combine contaminations and variants dictionary
allvariants = {**variants,**contaminations}

#the file called reads_name is opened and the length of the barcode-sequence is determined; the structure of the output-variable results is created: {sample:{barcode:variant=0}}
reads_input=""
print(my_dir)
objects=os.listdir(my_dir)
for filename in objects:
    if filename.endswith('gz'):
        reads=gzip.open(my_dir+filename, 'rb')
        reads_input=filename
        BCV_size=len(list(allvariants.keys())[0])
        read_count=0
        good_reads=0
        best_reads=0
        results={}
        print("\n"*10+"Sample being processed: %s" %reads_input)
        print("\n""Script is running...\n")
        for n in samples:
            results[samples[n]]={}
            for i in allvariants:
                results[samples[n]][allvariants[i]]=0
    else:
        print("no zipped file")
        continue
    print(results)
    
#for each line in the read_name-file, check if it is a DNA-sequence: If this is the case, store the sequence in line and increment read_count. If this is not the case, go to the next line
    print("checking for DNA sequences in file...sequences not matching ATGC are filtered out")
    for ln in reads:
        ln = str(ln).strip("b'\\n").upper()
        j+=1
        if j % 10000000 == 0:
            print(str(j/1000000) + " mio lines checked")
        if not ln.strip():
            continue
        for i in ln:
            if not i or i not in ["A","T","G","C"]:
                noDNA=True
                break
            noDNA=False
        if noDNA:
            continue
        line=ln.strip().upper()
        read_count+=1
#if the demultiplexing-variables are not defined in the config-file, BCS is always: ''
        if line[BCS1_size:BCS1_size+len(BCS1_right)]==BCS1_right.upper() and line[len(line)-BCS2_size-len(BCS2_left):len(line)-BCS2_size]==BCS2_left.upper():
                BCS=line[0:BCS1_size]+line[len(line)-BCS2_size:len(line)]
        else:
            continue

        #search for the the barcode-sequence of the read and store it in the variable BCV: If the constant hexamers are not found, the next read is analysed (no increment of good_reads or best_reads)
        # find is defined by arguments in brackets, result is the position, where the first argument (first in the bracket) is found
        #A is right constant hexamer, B is left constant hexamer. D and E are the reverse complement constant hexamers left and right.
        #rfind gives the position of the LAST occurence of the search term
        A=line.find(BCV_right,BCV_loc-BCV_margin+BCV_size)
        B=line.find(BCV_left,BCV_loc-BCV_margin-len(BCV_left))
        C=line.rfind(BCV_right,BCV_loc-BCV_margin+BCV_size,BCV_loc+BCV_margin+BCV_size+len(BCV_right))
        D=line.find(BCV_right_revcomp,BCV_loc_revcomp-BCV_margin+BCV_size)
        E=line.find(BCV_left_revcomp,BCV_loc_revcomp-BCV_margin-len(BCV_left))
        F=line.rfind(BCV_right_revcomp,BCV_loc_revcomp-BCV_margin+BCV_size,BCV_loc_revcomp+BCV_margin+BCV_size+len(BCV_right))
        if A-B==BCV_size+len(BCV_left):
            BCV=line[B+len(BCV_left):A]
        elif C-B==BCV_size+len(BCV_left):
            BCV=line[B+len(BCV_left):C]
        elif D-E==BCV_size+len(BCV_left):
            BCV=line[E+len(BCV_left):D]
        elif F-E==BCV_size+len(BCV_left):
            BCV=line[E+len(BCV_left):F] 
        else:
            continue
        
	#if the index sequence of the read can be matched with a sample stored in the variable samples, good_reads is incremented
        if BCS in samples:
            good_reads+=1
		#if the barcodesequence of the read can be matched with a variant stored in the variable variants, best_best_reads is incremented
            if BCV in allvariants:
                results[samples[BCS]][allvariants[BCV]]+=1
                best_reads+=1
       	#if the barcodesequence of the read is not find in the variant_file, the sequence is stored extra in results[samples] and incremented if this sequence is found repeatedly
            else:
                if BCV in results[samples[BCS]]:
                    results[samples[BCS]][BCV]+=1
                else:
                    results[samples[BCS]][BCV]=1

    output_file=my_dir+filename.replace('.gz','')
    #new_output_file=output_file.replace(".gz",".txt")
    f=open(output_file,'w')
    f.write("====== Generated with Barcode Script 1.4 ======\n\n")
    f.write("\nTotal number of reads: "+str(read_count))
    f.write("\nReads recovered: "+str(good_reads)+"\nReads with known variants: "+str(best_reads))
    f.write("\nReads with unknown variants: "+str(good_reads-best_reads)+"\n\nProportion of known variants:\n")
    print("\nTotal number of reads: "+str(read_count))
    print("Reads recovered: "+str(good_reads))
    print("Reads with known variants: "+str(best_reads))
    print("Reads with unknown variants: "+str(good_reads-best_reads)+"\n")
    for n in sorted(results,key=str.lower):
        if "output_file" in vars() and output_file!="":
            f.write("\nSample:"+n+" "*(16-len(n))+"# \n")

        for i in sorted(results[n],key=str.lower):
            if "output_file" in vars() and output_file!="":
                f.write(i+" "*(25-len(i)-len(str(results[n][i])))+str(results[n][i])+"\n")

    if "output_file" in vars() and output_file!="":
        f.write("\n")
        f.close()
    print("====== Script completed! ======")
