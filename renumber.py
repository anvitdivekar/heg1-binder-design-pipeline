import sys
from Bio import PDB
from Bio.PDB import PDBIO

in_path, out_path = sys.argv[1], sys.argv[2]
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure("design", in_path)
model = structure[0]

chain_b = model['B'] if 'B' in model else None
chain_a = model['A'] if 'A' in model else None
if chain_b is None or chain_a is None:
    raise SystemExit(f"Expected chains A and B, found: {[c.id for c in model]}")

binder_residues = [r for r in chain_b if PDB.is_aa(r)]
target_residues = [r for r in chain_a if PDB.is_aa(r)]

for i, res in enumerate(binder_residues, start=1):
    res.id = (' ', i, ' ')
n = len(binder_residues)
for i, res in enumerate(target_residues, start=n + 1):
    res.id = (' ', i, ' ')

io = PDBIO()
io.set_structure(structure)
io.save(out_path)
print(f"Wrote {out_path} — binder 1-{n}, target {n+1}-{n+len(target_residues)}")
