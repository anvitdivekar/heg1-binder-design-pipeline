import sys
from Bio import PDB

pdb_path, resfile_path = sys.argv[1], sys.argv[2]
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure("design", pdb_path)
model = structure[0]

with open(resfile_path, 'w') as f:
    f.write("ALLAA\nstart\n")
    for chain in model:
        if chain.id == 'A':
            for residue in chain:
                if PDB.is_aa(residue):
                    f.write(f"{residue.id[1]} A NATRO\n")
print(f"Wrote {resfile_path}")
