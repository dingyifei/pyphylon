#!/usr/bin/env python3
"""Build VFDB × phylon cross-reference table."""
import re
import pandas as pd

# 1. Parse VFDB FASTA headers
vfdb = {}
with open("data/external/VFDB/VFDB_setA_pro.fas") as f:
    for line in f:
        if line.startswith(">"):
            target_match = re.match(r">(\S+)", line)
            target_id = target_match.group(1).lstrip(">") if target_match else None

            vfg_match = re.match(r">(\w+)\(", line)
            vfg_id = vfg_match.group(1) if vfg_match else None

            gene_match = re.search(r"\)\s+\((\w+)\)", line)
            gene_name = gene_match.group(1) if gene_match else ""

            desc_match = re.search(r"\)\s+\(\w+\)\s+(.+?)\s+\[", line)
            description = desc_match.group(1) if desc_match else ""

            bracket_match = re.search(
                r"\[(.+?)\s+\(VF\d+\)\s*-\s*(.+?)\s+\(VFC\d+\)\]", line
            )
            vf_name = bracket_match.group(1) if bracket_match else ""
            vf_category = bracket_match.group(2) if bracket_match else ""

            species_match = re.search(r"\]\s+\[(.+?)\]\s*$", line)
            species = species_match.group(1) if species_match else ""

            vfdb[target_id] = {
                "vfg_id": vfg_id,
                "gene_name": gene_name,
                "description": description,
                "vf_name": vf_name,
                "vf_category": vf_category,
                "species": species,
            }

print(f"Parsed {len(vfdb)} VFDB entries")

# 2. Load BLAST results and phylon gene details
blast = pd.read_csv("output/data/5e_vfdb_blast_results.csv")
genes = pd.read_csv("output/data/5a_phylon_gene_details.csv")

print(f"BLAST results: {len(blast)} rows")
print(f"Gene details: {len(genes)} rows")

# 3. Strip allele suffix from query (CJejuni_C1010A0 -> CJejuni_C1010)
blast["gene"] = blast["query"].str.replace(r"A\d+$", "", regex=True)

# 4. Map target to VFDB annotations
blast["vf_name"] = blast["target"].map(
    lambda t: vfdb.get(t, {}).get("vf_name", "")
)
blast["vf_category"] = blast["target"].map(
    lambda t: vfdb.get(t, {}).get("vf_category", "")
)
blast["vfdb_gene"] = blast["target"].map(
    lambda t: vfdb.get(t, {}).get("gene_name", "")
)
blast["vfdb_description"] = blast["target"].map(
    lambda t: vfdb.get(t, {}).get("description", "")
)

print(f"\nVF categories found:")
for cat, count in blast["vf_category"].value_counts().items():
    print(f"  {cat}: {count}")

# 5. Join with phylon membership
merged = blast.merge(
    genes[["phylon", "gene", "is_exclusive", "product"]], on="gene", how="left"
)
print(f"\nMerged rows: {len(merged)}")
print(f"Matched to phylon: {merged['phylon'].notna().sum()}")
print(f"Missing phylon: {merged['phylon'].isna().sum()}")

if merged["phylon"].isna().any():
    unmatched = merged[merged["phylon"].isna()]["gene"].unique()[:10]
    print(f"Sample unmatched genes: {unmatched}")

# 6. Build phylon x virulence_category summary
pivot = (
    merged.dropna(subset=["phylon"])
    .groupby(["phylon", "vf_category"])
    .size()
    .unstack(fill_value=0)
)
print(f"\n=== Phylon x Virulence Category Table ===")
print(pivot.to_string())

# 7. Save outputs
merged.to_csv("output/data/vfdb_phylon_crossref.csv", index=False)
pivot.to_csv("output/data/vfdb_phylon_category_summary.csv")
print("\nSaved: output/data/vfdb_phylon_crossref.csv")
print("Saved: output/data/vfdb_phylon_category_summary.csv")

# 8. Per-phylon virulence details
print(f"\n=== Per-Phylon Virulence Details ===")
for phylon in sorted(merged["phylon"].dropna().unique()):
    subset = merged[merged["phylon"] == phylon]
    print(f"\n--- {phylon} ({len(subset)} VFDB hits) ---")
    for _, row in subset.iterrows():
        excl = "*" if row["is_exclusive"] else ""
        print(
            f"  {row['gene']}{excl}: {row['vfdb_gene']} - {row['vf_name']} "
            f"[{row['vf_category']}] ({row['identity']:.1f}%)"
        )
