from Bio.Data.IUPACData import protein_letters_3to1_extended as protein_letters_3to1, protein_letters_1to3_extended as protein_letters_1to3


def get_letters_3to1(x):
    '''
    convert 3 letter residue name x to 1 letter name
    '''
    x_new = x[0].upper() + x[1:].lower()
    if x_new in protein_letters_3to1:
        return protein_letters_3to1[x_new]
    else:
        return "X"
    

def get_letters_1to3(x):
    '''
    convert 1 letter residue name x to 3 letter name
    '''
    x_new = x.upper()
    return protein_letters_1to3[x_new]