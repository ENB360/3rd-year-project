data = ['a_pinguis_data', 'a_thaliana_data', 'b_oleracea_data', 'b_prasinos_data', 'b_vulgaris_data', 'c_annuum_data', 'c_arguta_data', 'c_crispus_data', 'c_lanatus_data', 'c_neogaea_data', 'c_papaya_data', 'c_paradoxa_data', 'c_reinhardtii_data', 'd_hygrometrica_data', 'e_guttata_data', 'e_siliculosus_data', 'g_max_data', 'h_brasiliensis_data', 'l_japonicus_data', 'm_domestica_data', 'm_paleacea_data', 'n_attenuata_data', 'n_nucifera_data', 'n_tabacum_data', 'o_rufipogon_data', 'p_dactylifera_data', 'p_glauca_data', 'p_patens_data', 'p_purpurea_data', 'r_communis_data', 's_bicolor_data', 't_aestivum_data', 't_timopheevii_data', 'v_faba_data', 'v_vinifera_data', 'z_mays_data']



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

# function for translation
def translate(seq):
    last_codon_start = len(seq) - 2
    protein = ""
    for start in range(0, last_codon_start,3):
        codon = seq[start:start+3]
        aa = Codon_table.get(codon, '')
        protein = protein + aa
    return(protein)


# function for hydrophobicity.
# WWW.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html col a.
# K&D Hydrophobicity values 
Hydrophobicity_table = {
    'A':'1.8', 'R':'-4.5', 'N':'-3.5', 'D':'-3.5', 'B':'-3.5', 'C':'2.5',
    'E':'-3.5', 'Q':'-3.5', 'Z':'-3.5', 'G':'-0.4', 'H':'-3.2', 'I':'4.5',
    'L':'3.8', 'K':'-3.9', 'M':'1.9', 'F':'2.8', 'P':'-1.6', 'S':'-0.8',
    'T':'-0.7', 'W':'-0.9', 'Y':'-1.3', 'V':'4.2', 'x':'0'
    }


# Function to get H_value
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


#Average Hydrophobicity value.
def avhval(proteinseq):
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
        ans1 = (v1 / v2)
        return (ans1)


# Function for length of amino acid sequence
def los(aa):
    full = len(aa)
    if aa[-1] == 'x':
        return(full - 1)
    else:
        return(full)
    


# Put it all together
file = open('filename')
f1 = file.read()
f2 = f1.split('\n')
Total_Hydrophobicities = []
check = [] 
GC_Percentages = []
Length_of = []
Average_Hydrophobicities = []
x = [f2[i:i+3] for i in range(0, len(f2),  3)]
for entry in x:
    name = entry[0]
    seq = entry[1]
    s2 = seq.split('sequence: ')[1]
    s3 = s2.lstrip()
    p = translate(s3)
    h = (h_val(p))
    h2 = (avhval(p))
    gc = (get_gc(s3))
    z = (los(p))
    print(name)
    print('Length of amino acid sequence = ' + str(z))
    print('Total amino acid hydrophobicity = ' + str(round(h, 2)))
    print('Average amino acid hydrophobicity = ' + str(round(h2, 2)))
    print('GC% = ' + str(round(gc, 2))+'%')
    print('\n')
    Length_of.append(float(z))
    Total_Hydrophobicities.append(float(h))
    Average_Hydrophobicities.append(float(h2))
    GC_Percentages.append(float(gc))
    n2 = name.split(':>')[0]
    check.append(n2)

print('All amino acid sequence lengths ' + str(Length_of))
print('\n')
print('All total hydrophobicities ' + str(Total_Hydrophobicities))
print('\n')
print('All average hydrophobicities ' + str(Average_Hydrophobicities))
print('\n')
print('All GC%\'s ' + str(GC_Percentages))
print('\n')
print('Gene lost in ' + str(set(data)^set(check)))



