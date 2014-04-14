#DNA_methods.py

def read_fasta_file(path):
    f= open(path,'r')
    seq_id = ((f.readline()).rstrip()).lstrip(">")
    sequence=""
    for line in f:
        sequence = sequence+line.rstrip()
    f.close()
    return (seq_id,sequence)

def reverse_complement(seq):
    seq=(seq[::-1]).upper()
    compliments = {'A':'T','C':'G','T':'A','G':'C'}
    revcom=""
    for i in seq:
        revcom=revcom+compliments[i]
    return revcom
    
def get_ORF(seq):
    seq=seq.upper()
    start=False
    stop=False
    rframe=0
    while (not(start) or not(stop)) and rframe<3: #loops until both a start and stop codon are found or until all 3 frames have been checked
        index = rframe
        start=False #resets start and stop for new frame
        stop=False        
        while (not(start) or not(stop)) and index<=len(seq)-3:#the inside loop, loops through 1 frame at a time until both codons are found or length of seq is exhuasted
            codon = seq[index:index+3]
            if codon == 'ATG':
                start=True
            if  start and (codon == 'TAG' or codon== 'TAA' or codon =='TGA'):#only finds a stop if a start has already been found
                stop=True
            index=index+3
        rframe+=1#moves onto next reading frame
      
    if start and stop:
        return rframe-1#rframe always one greater than it should be as it is incremented at end of loop
    else:
        return -1

def get_gene_by_ORF(seq,rf):
    seq=(seq[rf:]).upper()
    start= seq.index("ATG")
    stop_codons = ['TAG','TAA','TGA']
    gene=""
    codon=""
    while codon not in stop_codons:
        codon = seq[start:start+3]
        gene=gene+codon
        start+=3
    return gene

def translate(seq):
    amino_acids = {'ATG':'M','TGA':'stop','TAA':'stop','TAG':'stop','TTT':'F','TTC':'F','CTT':'L','CTC':'L','CTG':'L','CTA':'L','ATT':'I','ATC':'I','ATA':'I','GTA':'V','GTC':'V','GTG':'V','GTT':'V','TCA':'S','TCT':'S','TCG':'S','TCC':'S','CCA':'P','CCC':'P','CCT':'P','CCG':'P','ACT':'T','ACG':'T','ACC':'T','ACA':'T','GCG':'A','GCC':'A','GCA':'A','GCT':'A','TAT':'Y','TAC':'Y','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGA':'G','GGG':'G','GGC':'G','GGT':'G','CGA':'R','CGC':'R','CGG':'R','CGT':'R','TGT':'C','TGC':'C','TGG':'W','AGT':'S','AGC':'S'}
    seq=seq.upper()
    start=seq.index('ATG')
    codon=seq[start:start+3]
    protein=""
    while amino_acids[codon]!='stop'and start<=len(seq)-3:
        protein=protein+amino_acids[codon]
        start+=3
        codon=seq[start:start+3]
    protein=protein+amino_acids[codon]# adds on final stop outside of loop
    return protein
    
def get_fasta(seq_id,seq):
    fasta=seq_id+'\n'
    for i in range(len(seq)):
        fasta=fasta+seq[i]
        if (i+1)%60==0:
            fasta=fasta+'\n'
    return fasta
 

        
        
        
        
        
            
        
    
        
    
    
        
        
        

    
    
    