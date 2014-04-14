import sys
import DNA_methods

fasta_file=sys.argv[1]
protein_file=sys.argv[2]

seq_info=DNA_methods.read_fasta_file(fasta_file)

print "The reverse complement of the sequence is:\n"+DNA_methods.reverse_complement(seq_info[1])

ORF= DNA_methods.get_ORF(seq_info[1]) #returns ORF
gene= DNA_methods.get_gene_by_ORF(seq_info[1],ORF)# gets gene from the determined ORF
protein=DNA_methods.translate(gene)#translates the gene into amino acid sequence
fasta_protein=DNA_methods.get_fasta(seq_info[0]+' protein sequence',protein)

w=open(protein_file,'w')
w.write(fasta_protein)
w.close()
