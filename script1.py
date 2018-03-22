# Dictionary of codons and their corresponding amino acids. 
Codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'x', 'TAG':'x',
    'TGC':'C', 'TGT':'C', 'TGA':'x', 'TGG':'W',
    }

# Function for translation
def translate(seq):
    last_codon_start = len(seq) - 2
    protein = ""
    for start in range(0, last_codon_start,3):
        codon = seq[start:start+3]
        aa = Codon_table.get(codon, '')
        protein = protein + aa
    return(protein)




# Dictionary of hydrophobicity values. 
# WWW.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html
# Used the values found in column A.  
Hydrophobicity_table = {
    'A':'1.8', 'R':'-4.5', 'N':'-3.5', 'D':'-3.5', 'B':'-3.5', 'C':'2.5',
    'E':'-3.5', 'Q':'-3.5', 'Z':'-3.5', 'G':'-0.4', 'H':'-3.2', 'I':'4.5',
    'L':'3.8', 'K':'-3.9', 'M':'1.9', 'F':'2.8', 'P':'-1.6', 'S':'-0.8',
    'T':'-0.7', 'W':'-0.9', 'Y':'-1.3', 'V':'4.2',  'x':'0'
    }


    

# Function to get hydrophobicity value 
def h_val(proteinseq):
    full = len(proteinseq)
    d = []
    for start in range(full):
        proseq = proteinseq[start]
        val = Hydrophobicity_table.get(proseq, '')
        d.append(float(val))
    b = sum(d)
    return (b)


# Function to get GC content
def get_gc(dna):
    gcontent = dna.count('G')
    ccontent = dna.count('C')
    total = len(dna)
    if total == 0:
        return total
    else:
        ans = ((gcontent + ccontent) / total) * 100
        return (ans)




# Function to get average Hydrophobicity value.
def av_hval(proteinseq):
    full = len(proteinseq)
    f = []
    for start in range(full):
        proseq = proteinseq[start]
        val = Hydrophobicity_table.get(proseq, '')
        f.append(float(val))
    v1 = sum(f)
    v2 = len(f)
    if v2 == 0:
        return v2
    else:
        ans = (v1 / v2)
        return (ans)




# Put it all together
# Insert file name between '' on line below. 
file = open('filename') 
f1 = file.read()
f2 = f1.split('>lcl')
genes = []
for a in f2:
    a = '>' + a
    genes.append(a)
for gene in genes:
    indi = gene.split('\n')
    name = indi[0]
    seq = ''.join(indi[1:])
    print(name)
    print('DNA sequence: ' + seq)
    p = translate(seq)
    print('Amino acid sequence: ' + p)
    h = h_val(p)
    print('Total Hydrophobicity: '+ str(h))
    h2 = av_hval(p)
    print('Average Amino Acid Hydrophobicity: '+str(h2))
    GC = get_gc(seq)
    print ('GC Percentage: ' + str(GC)+'%')
    print('\n')

        



        
        
