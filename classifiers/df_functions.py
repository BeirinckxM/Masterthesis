def change_colnames(column):
    column = str(column)
    if 'Chain 1' in column:
        column = column.replace("Chain 1", "Heavy")
    if "Chain 2" in column:
        column = column.replace("Chain 2", "Light")
    else:
        column = column.replace("Description", "Epitope")
    return column

def change_organism_ID(data, organism_dict):
    heavy_org = []
    light_org = []
    for x in data['Heavy Species']:
        x = str(x)[:-2]
        try:
            x = organism_dict[x][0]
        except KeyError:
            x = x
        heavy_org.append(x)
    for x in data['Light Species']:
        x = str(x)[:-2]
        try:
            x = organism_dict[x][0]
        except KeyError:
            x = x
        light_org.append(x)
    data.insert(15, 'Heavy_Org', heavy_org)
    data.insert(45, 'Light_Org', light_org)

def get_Vgene(gene):
    gene = str(gene).replace('IGHV','')
    if len(gene) == 1:
        gene = '0' + gene
    else:
        pass
    
    return gene

def get_Vfamily(x):
    fam = str(x).split('*')[0]
    if len(fam) == 1:
        fam = '0' + fam
    
    return fam

def get_Jgene(gene):
    gene = str(gene).replace('IGHJ','')
    if len(gene) == 1:
        gene = '0' + gene
    else:
        pass
    
    return gene

def get_Jfamily(fam):
    if fam == 1:
        fam = '0' + fam
    else:
        pass
    
    return fam



def count_aromaticAA(seq):
    
    '''For a given sequence the aromatic amino acids (AAA) are counted. Function returns the total AAA
    count in the given sequence'''
    
    aromatic_AA = ['F', 'Y', 'W', 'H']
    counts = []
    for x in seq:
        count=0
        for i in x:
            if i in aromatic_AA:
                count +=1
        counts.append(count)
    return counts