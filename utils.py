import requests


def download_pdb(pdb_code, pdb_path, format = "cif"):
    '''
    download file from a pdb entry. supported format: pdb, cif, fasta
    '''
    if format == "fasta":
        r = requests.get(f"https://www.rcsb.org/fasta/entry/{pdb_code}")
    else:
        r = requests.get(f"https://files.rcsb.org/view/{pdb_code}.{format}")
    open(pdb_path, 'wb').write(r.content)
