from Bio import SeqIO
import os
import numpy
import sys


orderedHyperLoci = ["rl5a","rl6","rl12","rl13","ul1","ul9","ul11","ul20","ul73","ul74","ul120","ul139","ul146","us9"]
hypervar = set()
#Collect genotyoe and groups
genotypes = {}
infile = open("all_groups.txt")
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    hypervar.add(fields[0])
    if not (fields[0],fields[1]) in genotypes:
        genotypes[(fields[0],fields[1])] = fields[2]





#Collect cds
sequences = {}
geneScores = {}
unorderedGroups = {}
bestStrain = {}
outfile = open("genotypeClassification.txt","w")
outfile.write("Locus\tBestHit\n")
for seq_record in SeqIO.parse("cds.fasta","fasta"):
    locus = (str(seq_record.id)).lower()
    seq = str(seq_record.seq)

    if locus in hypervar:
        if not locus in unorderedGroups:
            unorderedGroups[locus] = ""
        geneScores[locus] = {}
        print "Analyzing hypervariable gene",locus
        tempFasta = open("toAnalyze.fasta","w")
        tempFasta.write(">"+locus+"\n"+seq+"\n")
        tempFasta.close()
        os.system("makeblastdb -dbtype nucl -in ./multiFasta/"+locus+".fa")
        os.system("blastn -query toAnalyze.fasta -db ./multiFasta/"+locus+".fa -outfmt 6 -out outputBlast.txt")
        print "Blast finished"
        infile = open("outputBlast.txt")
        bestScore = 0
        bestMatch = ""
        bestAlignmentLength = 0
        while True:
            line = infile.readline().rstrip()
            if not line:
                break
            fields = line.split("\t")
            if (len(seq) - (int(fields[7])-int(fields[6]))) < 50:
                if not fields[1] in bestStrain:
                    bestStrain[fields[1]] = 0
                bestStrain[fields[1]] += int(fields[11])



            if float(fields[11]) > bestScore:
                bestScore = float(fields[11])
                bestMatch = fields[1]
                bestAlignmentLength = int(fields[7])-int(fields[6])
        print bestAlignmentLength,len(seq)
        if len(seq) >2000:
            sys.stdin.read(1)
        if (len(seq) - bestAlignmentLength) > 50:
            print "Best match ",bestMatch,"is ",bestAlignmentLength," whereas the sequence is ",len(seq)
            sys.stdin.read(1)
        else:
            unorderedGroups[locus] = genotypes[(locus,bestMatch)]


for item in orderedHyperLoci:
    outfile.write(item+"\t"+unorderedGroups[item]+"\n")
bestOverallScore = 0
for item in bestStrain:
    if bestStrain[item] > bestOverallScore:
        bestOverallScore = bestStrain[item]
        bestHit = item

print "The best hit is ",bestHit,"with an overall score of ",bestOverallScore


            


            
            
            
            
            
            
            
            
            
            
            
            #if not fields[1] in geneScores[locus]:
            #    geneScores[locus][fields[1]] = float(fields[11])







#outfile = open("genotypeClassification.txt","w")
#outfile.write("Locus\t")
#for group in groups:
#    outfile.write(group+"\t")
#outfile.write("\n")
#for gene in geneScores:
#    print "Analyzing gene",gene
#    outfile.write(gene+"\t")
#    scoreRanking = {}
#    for group in groups:
#        numScores = 0
#        meanScore = 0
#        if not group in scoreRanking:
#            scoreRanking[group] = 0
#        for genotp in groups[group]:
#            if genotp in geneScores[gene]:
#                numScores +=1
#                meanScore += geneScores[gene][genotp]
#        meanScore = float(meanScore)/float(numScores)
#        #print "For genotype ",group,"the average is ",meanScore#
#        scoreRanking[group] = meanScore
#    scoreList = []
#    groupOrder = []
#    for group in scoreRanking:
#        groupOrder.append(group)
#        scoreList.append(scoreRanking[group])
#    newList = sorted(scoreList)
#    index = [newList.index(v) for v in scoreList]
#    print groupOrder
#    print scoreList
#    print index
#    for value in index:
#        outfile.write(str(value)+"\t")
#    outfile.write("\n")
#    sys.stdin.read(1)
    




    