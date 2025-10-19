# =============================================
# Viral Proteome Advanced Analysis Pipeline
# =============================================

import json
import os
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.SeqUtils import molecular_weight, IsoelectricPoint
from Bio.Seq import Seq
from sklearn.decomposition import PCA

# -----------------------------
# 1️⃣ File paths
# -----------------------------
json_file = "/mnt/d/Proteomics/proteome.json"
output_folder = "/mnt/d/Proteomics/Results"
os.makedirs(output_folder, exist_ok=True)

# -----------------------------
# 2️⃣ Load JSON
# -----------------------------
with open(json_file, "r", encoding="utf-8") as f:
    data = json.load(f)

proteins = data['results']
print(f"Total proteins loaded: {len(proteins)}")

# -----------------------------
# 3️⃣ Categories and counters
# -----------------------------
categories = ['envelope', 'membrane', 'capsid', 'other']
aa_counts_per_category = {cat: Counter() for cat in categories}

lengths, weights, pis, protein_categories = [], [], [], []

for protein in proteins:
    try:
        name = protein['proteinDescription']['recommendedName']['fullName']['value'].lower()
        seq_str = protein['sequence']['value']
        seq = Seq(seq_str)
    except KeyError:
        continue

    # Categorize
    matched = False
    for cat in categories[:-1]:
        if cat in name:
            aa_counts_per_category[cat].update(seq_str)
            protein_categories.append(cat)
            matched = True
            break
    if not matched:
        aa_counts_per_category['other'].update(seq_str)
        protein_categories.append('other')

    # Compute properties
    lengths.append(len(seq))
    weights.append(molecular_weight(seq, seq_type='protein'))
    ip_calc = IsoelectricPoint.IsoelectricPoint(str(seq))
    pis.append(ip_calc.pi())

# -----------------------------
# 4️⃣ Amino acid frequencies & enrichment
# -----------------------------
# Normalize
def normalize(counter):
    total = sum(counter.values())
    return {aa: count / total for aa, count in counter.items()}

aa_freq_per_category = {cat: normalize(cnt) for cat, cnt in aa_counts_per_category.items()}

# Total proteome
total_counts = Counter()
for seq in [p['sequence']['value'] for p in proteins]:
    total_counts.update(seq)
total_freq = {aa: count / sum(total_counts.values()) for aa, count in total_counts.items()}

# Enrichment
enrichment = {}
for cat, freq in aa_freq_per_category.items():
    enrichment[cat] = {aa: freq[aa] / total_freq.get(aa, 1e-6) for aa in freq}

df_enrichment = pd.DataFrame(enrichment).fillna(0)
df_freq = pd.DataFrame(aa_freq_per_category).fillna(0)

# Save CSVs
df_freq.to_csv(os.path.join(output_folder, "aa_frequencies_per_category.csv"))
df_enrichment.to_csv(os.path.join(output_folder, "aa_enrichment_per_category.csv"))
print("CSV files saved for frequencies and enrichment.")

# -----------------------------
# 5️⃣ Heatmap of enrichment
# -----------------------------
plt.figure(figsize=(10,6))
sns.heatmap(df_enrichment, annot=True, cmap="coolwarm", center=1)
plt.title("Amino Acid Enrichment per Category")
plt.xlabel("Category")
plt.ylabel("Amino Acid")
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "aa_enrichment_heatmap.png"))
plt.show()

# -----------------------------
# 6️⃣ 2D Proteome Mapping: MW vs pI
# -----------------------------
plt.figure(figsize=(10,6))
for cat in categories:
    cat_weights = [weights[i] for i,c in enumerate(protein_categories) if c==cat]
    cat_pis = [pis[i] for i,c in enumerate(protein_categories) if c==cat]
    plt.scatter(cat_pis, cat_weights, label=cat)
plt.xlabel("Isoelectric Point (pI)")
plt.ylabel("Molecular Weight (Da)")
plt.title("2D Proteome Map (MW vs pI)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "2D_proteome_map.png"))
plt.show()

# -----------------------------
# 7️⃣ PCA of amino acid frequencies
# -----------------------------
X = df_freq.T.values  # proteins as rows, aa as columns
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

plt.figure(figsize=(10,6))
plt.scatter(X_pca[:,0], X_pca[:,1])
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA of Amino Acid Frequencies")
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "PCA_aa_frequencies.png"))
plt.show()
print("PCA plot saved.")
