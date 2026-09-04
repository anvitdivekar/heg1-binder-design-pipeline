import glob, csv, sys, os

af2dir = sys.argv[1]
outpath = sys.argv[2]

rows = []
for scfile in glob.glob(os.path.join(af2dir, "*_scores.sc")):
    with open(scfile) as f:
        lines = [l for l in f if l.strip()]
    if len(lines) < 2:
        continue
    header = lines[0].split()[1:]
    for line in lines[1:]:
        vals = line.split()[1:]
        row = dict(zip(header, vals))
        rows.append(row)

if not rows:
    print("WARNING: no scored designs found — aggregate table is empty")
    sys.exit(0)

fieldnames = list(rows[0].keys())
with open(outpath, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(sorted(rows, key=lambda r: float(r.get('pae_interaction', 999))))

print(f"Wrote {outpath} — {len(rows)} scored designs")
