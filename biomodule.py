def reverseComplement(sequence):
    outSequence = ""
    accepted = ["A","C","T","G","a","c","t","g","N","n"]
    for a in range(len(sequence)-1,-1,-1):
        if sequence[a]=="A" or sequence[a]=="a":
            outSequence+="T"
        if sequence[a]=="T" or sequence[a]=="t":
            outSequence+="A"
        if sequence[a]=="G" or sequence[a]=="g":
            outSequence+="C"
        if sequence[a]=="C" or sequence[a]=="c":
            outSequence+="G"
        if sequence[a]=="N" or sequence[a]=="n":
            outSequence+="N"
        if not sequence[a] in accepted:
            print("Anomalous base found! Now exit....")
            exit()
    return outSequence