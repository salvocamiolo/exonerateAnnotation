from Bio import SeqIO
import biomodule as bm
import sys
import os

genomeName = sys.argv[1]
protName = sys.argv[2]

gffFile = open("annotation.gff","w")
warnFile = open("annotationWarnings.txt","w")
cdsFile = open("cds.fasta","w")
protFile = open("proteins.fasta","w")


for seq_record in SeqIO.parse(genomeName,"fasta"):
    genomeSeq = str(seq_record.seq)
    assemblyName = str(seq_record.id)


for seq_record in SeqIO.parse(protName,"fasta"):
    locus = str(seq_record.id)
    sequence = str(seq_record.seq)
    print "Analyzing locus",locus
    tempFasta = open("tempFasta.fasta","w")
    tempFasta.write(">"+locus+"\n"+sequence+"\n")
    tempFasta.close()
    os.system("exonerate --model protein2genome tempFasta.fasta "+genomeName+" --showtargetgff -s 0 -n 1  >outputExonerate")
    
    exResult = open("outputExonerate")
    line = exResult.readline().rstrip()
    while not "Query range:" in line:
        line = exResult.readline().rstrip()
        if line is None:
            print "WARNING no protein found for ", locus
            warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
 
    split1 = line.split(": ")
    split2 = split1[1].split(" -> ")


    exResult.close()
    print float(int(split2[1])-int(split2[0]))/float(len(sequence)*3)
    if float(int(split2[1])-int(split2[0]))/float(len(sequence)) < 0.9 or int(split2[0])>0:
        print "Applying refinement"
        os.system("exonerate --model protein2genome tempFasta.fasta "+genomeName+" --showtargetgff -s 0 -n 1 --refine full >outputExonerate")













    exResult = open("outputExonerate")
    line = exResult.readline().rstrip()
    while not line == "# --- START OF GFF DUMP ---":
        line = exResult.readline().rstrip()
        if line is None:
            #print "WARNING no protein found for ", locus
            warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
    for a in range(10):
        line = exResult.readline().rstrip()

    gene = {}
    exon = {}
    warnings = []
    #Collect the exonerate output
    while not line == "# --- END OF GFF DUMP ---":
        line = exResult.readline().rstrip()
        if line is None:
            #print "WARNING no protein found for ", locus
            warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
        fields = line.split("\t")
        if not line == "# --- END OF GFF DUMP ---":
            if fields[2]=="gene":
                if not locus in gene:
                    gene[locus] = (fields[3],fields[4],fields[6])
                else:
                    if (int(fields[4]) - int(fields[3])) > (int(gene[locus][1]) - int(gene[locus][0])):
                        gene[locus] = (fields[3],fields[4],fields[6])
            
            if fields[2] == "exon":
                if not locus in exon:
                    exon[locus] = []
                exon[locus].append((fields[3],fields[4],fields[6]))
                if "frameshifts" in fields[8]:
                    warnings.append("Frameshifts in exon "+str(fields[3])+" "+str(fields[4])+" "+fields[8])
   
    newList = sorted(exon[locus], key=lambda x: x[1])
    exon[locus] = newList
    

    
    #Check the coding sequence
    cdsSeq = ""
    if exon[locus][0][2]=="+":
        for item in exon[locus]:
            cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
        #print "CDS sequence"
        #print cdsSeq
        #print exon[locus]

        #Check start codon
        if not cdsSeq[:3]=="ATG":
            foundCodon = 0
            for a in range(len(sequence)-len(cdsSeq)/3+100):
                #print "New start codons"
                #print genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3]
                if genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3] == "ATG":
                    exon[locus][0]=(int(exon[locus][0][0])-a*3-3,exon[locus][0][1],exon[locus][0][2])
                    gene[locus] = (int(exon[locus][0][0])-a*3-3,gene[locus][1],gene[locus][2])
                    cdsSeq = ""
                    for item in exon[locus]:
                        cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]     
                    #print "New modified cds"
                    #print exon[locus]     
                    #print cdsSeq
                    foundCodon = 1
                    break
            if foundCodon == 1:
                warnFile.write("WARNING A valid start codon was found "+str(a*3)+" codons downstream the predicted gene "+locus+"\n")
            else:
                warnFile.write("CRITICAL! No Start codon for gene "+locus+"\n")

        

            

            


        #Check stop codon
        if not cdsSeq[-3:]=="TGA" and not cdsSeq[-3:]=="TAA" and not cdsSeq[-3:]=="TAG":
            foundCodon = 0
            for a in range(len(sequence)-len(cdsSeq)/3+100):
                #print "New Stop codons"
                #print genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3]
                newStop = genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3]
                if newStop == "TAA" or newStop=="TGA" or newStop=="TAG":
                    exon[locus][0]=(int(exon[locus][0][0]),int(exon[locus][0][1])+a*3+3,exon[locus][0][2])
                    gene[locus] = (gene[locus][0], int(exon[locus][0][1])+a*3+3, gene[locus][2])
                    cdsSeq = ""
                    for item in exon[locus]:
                        cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]     
                    #print "New modified cds"
                    #print cdsSeq
                    #print exon[locus]
                    foundCodon = 1
                    break
            if foundCodon == 1:
                warnFile.write("WARNING A valid stop codon was found "+str(a*3)+" codons upstream the predicted gene "+locus+"\n")
            else:
                warnFile.write("CRITICAL! No Stop codon for gene "+locus+"\n")

        gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(gene[locus][0])+"\t"+str(gene[locus][1])+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
        gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(gene[locus][0])+"\t"+str(gene[locus][1])+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
        numExon = 1
        for item in exon[locus]:
            gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
            numExon += 1
        #sys.stdin.read(1)
        cdsFile.write(">"+locus+" +\n"+cdsSeq+"\n")

    else:
        for a in range(len(exon[locus])-1,-1,-1):
            #print "Starting pos",int(exon[locus][a][0])-1
            #print genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ]
            cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ])
        #print "Final CDS"
        #print cdsSeq
        #print exon[locus]
        #print locus
        #sys.stdin.read(1)

     #Check start codon
        if not cdsSeq[:3]=="ATG":
            foundCodon = 0
            for a in range(len(sequence)-len(cdsSeq)/3+100):
                #print "New start codons"
                newStart = bm.reverseComplement(genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3])
                #print newStart
                if newStart == "ATG":
                    exon[locus][-1]=(int(exon[locus][-1][0]),int(exon[locus][-1][1])+a*3+3,exon[locus][-1][2])
                    gene[locus] = (gene[locus][0],int(exon[locus][-1][1])+a*3+3,gene[locus][2])
                    cdsSeq = ""
                    for a in range(len(exon[locus])-1,-1,-1):
                        #print "Starting pos",int(exon[locus][a][0])-1
                        cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ])   
                    #print "New modified cds"
                    #print exon[locus]     
                    #print cdsSeq
                    #sys.stdin.read(1)
                    foundCodon = 1
                    break
            if foundCodon == 1:
                warnFile.write("WARNING A valid start codon was found "+str(a*3)+" codons upstream the predicted start site\n")
            else:
                warnFile.write("CRITICAL! No Start codon for gene "+locus+"\n")

        #Check stop codon
        if not cdsSeq[-3:]=="TGA" and not cdsSeq[-3:]=="TAA" and not cdsSeq[-3:]=="TAG":
            foundCodon = 0
            for a in range(len(sequence)-len(cdsSeq)/3+100):
                #print "New Stop codons"
                newStop = bm.reverseComplement(genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3])
                #print newStop
                if newStop == "TAA" or newStop=="TGA" or newStop=="TAG":
                    exon[locus][0]=(int(exon[locus][0][0])-a*3-3,exon[locus][0][1],exon[locus][0][2])
                    gene[locus] =(int(exon[locus][0][0])-a*3-3, gene[locus][1],gene[locus][2])
                    cdsSeq = ""
                    for b in range(len(exon[locus])-1,-1,-1):
                        #print "Starting pos",int(exon[locus][a][0])-1
                        cdsSeq+=bm.reverseComplement(genomeSeq[ int(exon[locus][b][0])-1: int(exon[locus][b][1]) ])    
                    #print "New modified cds"
                    #print cdsSeq
                    #print exon[locus]  
                    #print a
                    #sys.stdin.read(1)   
                    foundCodon = 1
                    break
            if foundCodon == 1:
                warnFile.write("WARNING A valid stop codon was found "+str(a*3)+" codons upstream the predicted gene "+locus+"\n")
            else:
                warnFile.write("CRITICAL! No Stop codon for gene "+locus+"\n")

        gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(gene[locus][0])+"\t"+str(gene[locus][1])+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
        gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(gene[locus][0])+"\t"+str(gene[locus][1])+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
        numExon = 1
        for item in exon[locus]:
            gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
            numExon += 1
        
        #sys.stdin.read(1)

        cdsFile.write(">"+locus+" -\n"+cdsSeq+"\n")


gffFile.close()
for item in warnings:
    warnFile.write(item+"\n")

warnFile.close()








    





    

    
    
    
    
    

    