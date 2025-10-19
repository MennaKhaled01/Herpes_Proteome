# Herpes_Proteome
This repository presents an advanced bioinformatics analysis of the Human herpesvirus 1 (strain 17) proteome (UniProt ID: UP000009294), quantifying how protein localization (Envelope, Membrane, Capsid) correlates with molecular signatures such as amino acid composition, molecular weight (MW), and isoelectric point (pI)

# Pipeline Overview

The analysis pipeline systematically processes the HSV-1 proteome to evaluate molecular and compositional trends across functional categories.
The workflow includes the following stages:

1-Protein Categorization: Each protein is assigned to a structural category based on annotation keywords.

2-Physicochemical Property Calculation: Molecular Weight (MW) and Isoelectric Point (pI) are computed for every protein.

3-Amino Acid Analysis: Frequencies and enrichment ratios (category-specific vs. global proteome average) are calculated.

4-Data Visualization:

2D Proteome Map — MW vs. pI, colored by category.

Amino Acid Enrichment Heatmap — compositional deviations across categories.

PCA Plot — clustering of protein categories based on AA composition.

# Quantitative Results and Key Findings

**A. Amino Acid Enrichment Analysis: Structural Signatures**
The enrichment analysis identifies significant compositional biases required by structural constraints, with the $\text{Membrane}$ category showing the most extreme deviation from the proteome average.
​The amino acid enrichment analysis revealed clear compositional distinctions among the HSV-1 protein categories, with particularly striking patterns observed in membrane-associated proteins. Tryptophan (W) exhibited the highest enrichment in membrane proteins, approximately 89.2% above the proteome average, consistent with its well-established role in anchoring proteins to lipid bilayers. Similarly, Glycine (G) was strongly enriched in the membrane category—around 67.6% above average, indicating its contribution to flexibility within transmembrane helices. In contrast, Asparagine (N) showed a marked depletion in membrane proteins (~58.3% below average) while being enriched in capsid proteins (~32.2% above average), reflecting divergent structural requirements between hydrophobic and soluble environments. Finally, Lysine (K) was the most depleted residue in membrane proteins (~59.4% below average), aligning with the exclusion of positively charged residues from the hydrophobic core of lipid membranes. Collectively, these findings underscore the distinct biochemical constraints shaping protein composition across different viral compartments.
<img width="1000" height="600" alt="aa_enrichment_heatmap" src="https://github.com/user-attachments/assets/439618e2-2c4c-4b52-99ad-fcb0f37ee10e" />

**B. Distribution of Amino Acids (Absolute Frequencies)**

Absolute amino acid frequency analysis revealed distinct compositional patterns across the four structural categories of HSV-1 proteins. Alanine (A) displayed the highest overall abundance in membrane proteins (15.35%), exceeding its frequency in all other categories and highlighting its structural role in forming compact hydrophobic regions. Similarly, Glycine (G) showed a strong enrichment in membrane proteins (13.19%) compared to capsid (8.20%) and envelope (7.58%) proteins, reflecting its importance in providing flexibility to transmembrane helices. Leucine (L), another hydrophobic residue, was notably more frequent in membrane proteins (12.23%) than in other categories, reinforcing the dominance of nonpolar residues in lipid-interacting domains. Meanwhile, Proline (P) exhibited moderate enrichment in membrane proteins (10.31%) and envelope proteins (9.26%), suggesting its contribution to helix bending and structural diversity. Collectively, these data confirm that HSV-1 membrane proteins are overwhelmingly composed of small, nonpolar amino acids—an adaptive signature of their integration within the viral envelope’s hydrophobic environment.
<img width="1000" height="600" alt="aa_heatmap" src="https://github.com/user-attachments/assets/57cf1bc1-6655-4c80-abd3-03bd85004c60" />

The combined proportion of non-polar/small amino acids (A, G, L, P) in Membrane proteins reaches 51.08%, compared to 36.69% in Envelope proteins. This highlights the strong hydrophobic bias imposed by the lipid bilayer, which shapes membrane protein evolution and topology.

**C. Biophysical and Compositional Separation**
The 2D proteome map illustrates the distribution of HSV-1 proteins based on their molecular weight (MW) and isoelectric point (pI), offering a comprehensive view of proteomic diversity. Molecular weights range from approximately 10 kDa to over 300 kDa, reflecting the wide structural variability among viral proteins. The pI values span from 4.5 to 10, indicating substantial diversity in protein charge characteristics. Notably, capsid proteins predominantly cluster at higher molecular weights and moderate pI values**, while membrane proteins tend to be smaller and more alkaline, consistent with their compact, hydrophobic nature and functional localization within the viral envelope.

<img width="1000" height="600" alt="2D_proteome_map" src="https://github.com/user-attachments/assets/71cc532a-7a69-4580-9c88-0b1e963ee6b1" />

# PCA of Amino Acid Frequencies
Principal Component Analysis (PCA) was performed to determine the major compositional axes differentiating HSV-1 protein categories. The first two components captured a substantial proportion of the variance, with PC1 explaining approximately 47% and PC2 accounting for around 18% of the total variation. Membrane proteins exhibited the most pronounced separation along PC1, primarily driven by the contrast between hydrophobic and hydrophilic amino acid content. In contrast, the Envelope, Capsid, and Other protein groups clustered more closely together, indicating a higher degree of compositional similarity among the non-membranous proteins.The PCA confirms that hydrophobic enrichment is the dominant compositional feature differentiating HSV-1 protein classes, underpinning their distinct structural and environmental adaptations.

<img width="1000" height="600" alt="PCA_aa_frequencies" src="https://github.com/user-attachments/assets/d10a6df9-0e52-42f6-a413-ddbe65efc629" />
