# KRAS G12C Bioinformatics Pipeline

*Sequence Retrieval · BLAST Analysis · Structural Biology · AlphaFold Validation · Molecular Docking*

---


## Colab Notebook Link
[Open in Google Colab](https://colab.research.google.com/drive/1BP86eW18Xxtl44SkRa-XR9PsJX-8MTrZ?usp=sharing)

## 1. Project Identity

| Field | Details |
|---|---|
| **Target gene** | KRAS (*Kirsten Rat Sarcoma Viral Proto-oncogene*) — *Homo sapiens* |
| **Mutation focus** | G12C oncogenic point mutation |
| **Drug of interest** | Sotorasib / AMG 510 — PubChem CID: 137278711 |
| **Experimental structure** | PDB 6OIM (KRAS G12C–Sotorasib co-crystal, X-ray crystallography) |
| **Predicted structure** | AlphaFold AF-P01116-F1 (EBI AlphaFold Database, model v6) |
| **Docking engine** | GNINA v1.1 (CNN-augmented AutoDock Vina, built Dec 18 2023) |
| **Runtime environment** | Google Colaboratory — Python 3, Ubuntu 22.04 |
| **Reproducibility seed** | `--seed 42` (fixed across all docking runs) |

---

## 2. Background

KRAS (*Kirsten Rat Sarcoma Viral Proto-oncogene*) is one of the most frequently mutated genes in human cancer, implicated in non-small cell lung cancer (NSCLC), colorectal cancer, and pancreatic ductal adenocarcinoma. The **G12C point mutation** — a single nucleotide substitution replacing glycine with cysteine at position 12 — locks the protein in a constitutively active GTP-bound state, driving uncontrolled cell proliferation through downstream RAS–MAPK and PI3K–AKT signalling pathways.

**Sotorasib (AMG 510, brand name Lumakras)** is the first FDA-approved covalent inhibitor targeting KRAS G12C. Its acrylamide warhead forms a selective, irreversible thioether bond with the mutant cysteine 12 within the **switch II pocket (S-IIP)** — a cryptic site accessible only in the GDP-bound inactive conformation. This covalent, mutation-selective mechanism makes Sotorasib the ideal candidate for this structure-based computational study.

---

## 3. Pipeline Overview

The pipeline consists of nine sequential Google Colab cells, each handling one distinct biological or computational task:

```
Cell 0   GNINA binary download & validation
   │
Cell 1   Environment setup & project folder creation
   │
Cell 2   Project IDs, accession numbers & source links
   │
Cell 3   Sequence retrieval → transcription → translation → hotspot check
   │
Cell 4   BLASTn homology search (NCBI nt database)
   │
Cell 5   Structure retrieval (PDB 6OIM + AlphaFold) & receptor cleaning
   │
Cell 6   Structural validation — RMSD superimposition + pLDDT analysis
   │
Cell 7   Ligand preparation — Sotorasib SDF download & format conversion
   │
Cell 8   GNINA molecular docking (10 poses, exhaustiveness 8, seed 42)
   │
Cell 9   Results extraction, scoring table, plots & 3D visualisation
```
# Environment Setup

| Tool / Library | Type | Purpose |
|---|---|---|
| Google Colab | Platform | Cloud-based Jupyter notebook environment with Python 3.10 and GPU (T4) support |
| GNINA v1.1 | System tool | Molecular docking engine that scores binding poses using physics-based energy and a convolutional neural network |
| Open Babel | System tool | Converts Sotorasib between chemical file formats (SDF → PDB → PDBQT) and optimises 3D atomic coordinates |
| wget | System tool | Downloads the GNINA binary from GitHub release links during setup |
| zip | System tool | Compresses the entire project output into a single downloadable ZIP archive |
| biopython | Python library | Fetches KRAS sequence from NCBI, parses GenBank/FASTA/PDB files, runs BLAST, and calculates RMSD between structures |
| requests | Python library | Downloads files from external databases including NCBI, UniProt, AlphaFold, PubChem, and RCSB PDB |
| pandas | Python library | Organises results into structured tables and saves them as CSV files |
| matplotlib | Python library | Generates and saves the pLDDT confidence plot and the docking affinity bar chart |
| py3Dmol | Python library | Renders interactive 3D visualisations of protein structures and the docked ligand inside Colab |
| rdkit | Python library | Chemistry toolkit included for molecular structure support |
| os | Standard library | Creates project folders, checks file existence, and manages file paths |
| re | Standard library | Extracts docking score values per pose from the GNINA compressed SDF output |
| gzip | Standard library | Reads the compressed GNINA docking output file (.sdf.gz) without manual decompression |
| subprocess | Standard library | Runs the GNINA binary from within Python and captures output to verify it works correctly |
| shutil | Standard library | Creates the final ZIP archive of the complete project folder for download |

---

## 4. Methodology

### 4.1 Sequence Retrieval & Processing

**Cells 2–3 | Tools: Biopython Entrez, SeqIO | Databases: NCBI Nucleotide, NCBI Protein, UniProt**

The pipeline begins by programmatically fetching the canonical human KRAS gene record from NCBI. The **GenBank format** is specifically requested — rather than plain FASTA — because GenBank records contain structured feature annotations including the CDS (coding sequence) region with precise nucleotide start and stop positions.

The CDS feature is located by iterating through all annotated features and selecting the entry where `feature.type == "CDS"`. The extracted coding DNA is then processed through two sequential steps reflecting the **central dogma of molecular biology**:

```
Step 1 — Transcription:   Coding DNA (567 bp)  →  mRNA (567 nt)      via .transcribe()
Step 2 — Translation:     mRNA (567 nt)         →  Protein (188 aa)   via .translate(to_stop=True)
```

**Hotspot verification**

After translation, three clinically critical oncogenic hotspot positions are verified by direct index lookup:

| Hotspot | Array index | Wild-type residue | Known oncogenic variants |
|---|---|---|---|
| G12 | `protein[11]` | Glycine (G) | G12C · G12D · G12V · G12S · G12R · G12A |
| G13 | `protein[12]` | Glycine (G) | G13D · G13C |
| Q61 | `protein[60]` | Glutamine (Q) | Q61H · Q61L · Q61R · Q61K |

Two additional protein sequences are fetched in parallel as independent references: **NP_004976.2** from NCBI Protein and **P01116** from the UniProt REST API.

**Confirmed output**

| Sequence | Source | Length |
|---|---|---|
| Coding DNA | NM_004985.5 | 567 bp |
| mRNA | Transcribed from coding DNA | 567 nt |
| Translated protein | Derived from mRNA | 188 aa |
| Reference protein | NP_004976.2 (NCBI Protein) | 188 aa |
| Canonical protein | P01116 (UniProt) | 188 aa |
| Hotspot G12 | Confirmed wild-type | Glycine (G) |
| Hotspot G13 | Confirmed wild-type | Glycine (G) |
| Hotspot Q61 | Confirmed wild-type | Glutamine (Q) |

---

### 4.2 BLAST Homology Search

**Cell 4 | Tool: NCBIWWW.qblast, NCBIXML | Database: NCBI nt (non-redundant nucleotide)**

A **nucleotide BLAST (BLASTn)** search is executed against the NCBI `nt` database using the 567 bp KRAS coding DNA as the query. The search is submitted via `NCBIWWW.qblast` and results are returned in XML format. The top 5 hits are requested using `hitlist_size=5`.

For each alignment, only the **best HSP (High-Scoring Segment Pair)** is analysed. Percentage identity is computed as:

```
identity (%) = (hsp.identities / hsp.align_length) × 100
```

Results are saved in three formats: raw XML for archival, a parsed CSV, and a PNG visual summary table.

**Results**

| Rank | Hit | Identity (%) | E-value | Alignment length |
|---|---|---|---|---|
| 1 | *Homo sapiens* — PV477844.1 | 100.0 | 0.0 | 567 |
| 2 | *Pan troglodytes* — XM_054664451.2 | 100.0 | 0.0 | 567 |
| 3 | *Pan paniscus* — XM_054443149.2 | 100.0 | 0.0 | 567 |
| 4 | *Pan troglodytes* — XM_528758.9 | 100.0 | 0.0 | 567 |
| 5 | *Pan paniscus* — XM_054443148.2 | 100.0 | 0.0 | 567 |

All five hits returned **100% identity at E-value 0.0** across the full 567 bp alignment — the best possible BLAST result. Hits 2–5 are chimpanzee and bonobo sequences, reflecting the extreme evolutionary conservation of KRAS across primates.

> **Runtime note:** `NCBIWWW.qblast` submits to remote NCBI servers and may take 5–15 minutes depending on queue load. For high-throughput use, the standalone BLAST+ toolkit with a local `nt` database is recommended.

---

### 4.3 Structure Retrieval & Receptor Preparation

**Cell 5 | Tools: requests, Biopython PDBParser + PDBIO | Sources: RCSB PDB, AlphaFold EBI**

Two PDB structure files are obtained, serving complementary roles:

**Experimental structure — PDB 6OIM**
Downloaded from RCSB PDB (`https://files.rcsb.org/download/6OIM.pdb`). This is an X-ray crystal structure of KRAS G12C in complex with Sotorasib, providing the experimentally resolved binding pose within the switch II pocket. It serves as the docking receptor and the reference for binding box centroid calculation.

**Predicted structure — AlphaFold AF-P01116-F1**
Downloaded from the EBI AlphaFold Protein Structure Database. Model v6 is attempted first, with a deterministic fallback chain (v6 → v5 → v4 → v3 → RCSB-hosted). This represents the wild-type apo KRAS conformation, used for structural comparison and confidence analysis.

**Non-protein content of PDB 6OIM**

| Residue code | Description | Atom count |
|---|---|---|
| HOH | Crystallographic water molecules | 207 |
| MOV | Sotorasib / AMG 510 (co-crystallised ligand) | 41 |
| GDP | Guanosine diphosphate (endogenous nucleotide) | 28 |
| MG | Magnesium ion (cofactor) | 1 |

**Receptor cleaning**

All HETATM records are stripped using a custom Biopython `Select` subclass, retaining only standard amino acid residues:

```python
class KeepProteinOnly(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "   # blank insertion code = standard amino acid
```

This cleaned, protein-only PDB file is the direct receptor input for GNINA docking.

---

### 4.4 Structural Validation: RMSD & pLDDT Analysis

**Cell 6 | Tools: Biopython Superimposer, matplotlib, py3Dmol**

Both structures are first visualised interactively in 3D using py3Dmol — the experimental 6OIM as a spectrum-coloured cartoon with heteroatoms as green sticks, and the AlphaFold model separately for direct visual comparison. A quantitative structural analysis is then performed through two methods.

**RMSD superimposition**

Cα atoms are extracted from both structures by iterating through all models, chains, and residues, retaining entries containing a `CA` atom key. The shorter list determines the count of aligned atoms for a fair pairwise comparison:

```python
n = min(len(af_ca), len(exp_ca))
sup = Superimposer()
sup.set_atoms(exp_ca[:n], af_ca[:n])   # aligns AlphaFold onto experimental
rmsd = sup.rms
```

**pLDDT confidence extraction**

AlphaFold stores per-residue confidence (pLDDT) in the B-factor column of its PDB output. This is read directly from each Cα atom and plotted as a per-residue line chart with three reference thresholds:

| pLDDT threshold | Confidence level | Meaning |
|---|---|---|
| 90 | Very high | Backbone expected to be highly accurate |
| 70 | Confident | Generally correct overall structure |
| 50 | Low | Potentially disordered or flexible region |

**Results**

| Metric | Value | Interpretation |
|---|---|---|
| Aligned Cα residues | 167 | Sufficient for meaningful comparison |
| RMSD (AlphaFold vs 6OIM) | 4.807 Å | Moderate divergence — expected due to ligand-induced conformation |
| Average pLDDT | 91.52 / 100 | Very high confidence — backbone highly reliable |

The 4.807 Å RMSD reflects the structural difference between the ligand-bound switch II pocket (6OIM) and the apo state (AlphaFold) — a well-documented conformational change central to the KRAS G12C mechanism. The average pLDDT of 91.52 confirms the AlphaFold model is a reliable structural reference.
<img src="C:\Users\Rameen_Sajjad\Downloads\KRAS_Bioinformatics_Pipeline\figures" alt="AlphaFold PLDDT Plot" width="500"/>

---

### 4.5 Ligand Preparation

**Cell 7 | Tools: requests, Open Babel (obabel) | Source: PubChem REST API**

Sotorasib is retrieved from the PubChem Compound REST API by CID (137278711). The 3D SDF endpoint is requested first; a standard 2D SDF is used as fallback if the 3D version is unavailable:

```
Primary:  https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/137278711/SDF?record_type=3d
Fallback: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/137278711/SDF
```

Open Babel converts the SDF to two required formats using the `--gen3d` flag to ensure valid 3D geometry:

```bash
# For 3D visualisation in py3Dmol
obabel sotorasib_pubchem_3d.sdf -O sotorasib.pdb --gen3d

# For GNINA docking input (AutoDock PDBQT format)
obabel sotorasib_pubchem_3d.sdf -O sotorasib.pdbqt --gen3d
```

> **Open Babel kekulization warning:** Sotorasib's complex fused heterocyclic ring system may trigger a non-fatal warning: *"Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders"*. This does not affect the 3D coordinates or atom connectivity used in docking.

---

### 4.6 GNINA Molecular Docking

**Cell 8 | Tool: GNINA v1.1 | Method: CNN-augmented AutoDock Vina | 10 poses, exhaustiveness 8**

GNINA extends AutoDock Vina's force-field scoring with a convolutional neural network (CNN) layer trained on experimental protein–ligand crystal structures. This produces two independent scoring signals per pose: a physics-based binding estimate (`minimizedAffinity`) and a learned pose quality score (`CNNscore`).

**Docking box definition**

The search box is derived automatically from the co-crystallised ligand (MOV) coordinates in 6OIM. All MOV HETATM atoms are extracted and their centroid is computed as the arithmetic mean of Cartesian coordinates:

```python
center_x = sum(c[0] for c in coords) / len(coords)
center_y = sum(c[1] for c in coords) / len(coords)
center_z = sum(c[2] for c in coords) / len(coords)
box_size  = 22  # Å, uniform in X, Y, and Z
```

**GNINA command**

```bash
./gnina \
  -r kras_6OIM_clean_receptor.pdb \
  -l sotorasib_pubchem_3d.sdf \
  --center_x <X>  --center_y <Y>  --center_z <Z> \
  --size_x 22     --size_y 22     --size_z 22 \
  --exhaustiveness 8 \
  --num_modes 10 \
  --seed 42 \
  -o kras_sotorasib_gnina_docked.sdf.gz
```

**Docking parameters**

| Parameter | Value | Rationale |
|---|---|---|
| `--exhaustiveness` | 8 | Standard search depth; increase to 16–32 for publication-grade precision |
| `--num_modes` | 10 | Number of distinct binding poses returned |
| `--seed` | 42 | Fixed random seed for full reproducibility |
| Box size | 22 Å × 22 Å × 22 Å | Covers the switch II pocket with adequate sampling margin |
| Receptor | Protein-only clean PDB | All HETATM atoms removed before docking |
| Ligand | 3D SDF from PubChem | GNINA reads SDF directly |
| Output | `.sdf.gz` | All 10 poses in a single compressed SDF file |

**Score fields extracted**

| Field | Type | Description |
|---|---|---|
| `minimizedAffinity` | Force-field (kcal/mol) | Primary docking score; more negative = stronger predicted binding |
| `CNNscore` | Neural network (0–1) | Pose quality; higher = more physically realistic binding mode |
| `CNNaffinity` | Neural network | CNN-predicted affinity, complementary to force-field score |
| `CNNaffinity_variance` | Neural network | Uncertainty estimate for CNN affinity |

**Docking results — all 10 poses**

| Pose | Affinity (kcal/mol) | Intramol. (kcal/mol) | CNN pose score | CNN affinity |
|---|---|---|---|---|
| **1 ★** | **−10.67** | 1.08 | **0.3321** | **7.251** |
| 2 | −10.01 | 2.76 | 0.2542 | 7.072 |
| 3 | −6.10 | 0.19 | 0.2269 | 6.483 |
| 4 | −6.20 | 1.46 | 0.1916 | 6.510 |
| 5 | −6.48 | −0.29 | 0.1896 | 6.569 |
| 6 | −6.78 | −1.15 | 0.1434 | 6.804 |
| 7 | −6.70 | 0.44 | 0.1415 | 6.567 |
| 8 | −9.66 | 1.04 | 0.1329 | 6.896 |
| 9 | −7.57 | 0.78 | 0.1324 | 6.546 |
| 10 | −8.10 | −1.71 | 0.1302 | 6.555 |

★ Best pose — selected by minimum `minimizedAffinity`

**Affinity benchmarks**

| Affinity range (kcal/mol) | Predicted binding strength |
|---|---|
| > −6 | Weak or non-specific |
| −6 to −8 | Moderate |
| −8 to −10 | Strong |
| < −10 | Very strong — **Pose 1: −10.67 kcal/mol** |

Pose 1 achieves both the most negative affinity and the highest CNNscore across all 10 poses. This convergence of force-field and neural network evidence makes it the most confident predicted binding mode, consistent with Sotorasib's documented sub-nanomolar IC₅₀ against KRAS G12C.

> **GPU note:** Without a GPU runtime, GNINA warns *"No GPU detected. CNN scoring will be slow."* and falls back to CPU-only mode. All results remain correct — only runtime is affected (3–8 min GPU vs 15–30 min CPU).

---

### 4.7 Results Extraction & Analysis

**Cell 9 | Tools: gzip, re, pandas, matplotlib, py3Dmol, Open Babel**

The compressed GNINA output (`.sdf.gz`) is decompressed in memory using Python's `gzip` module. Score values are extracted with regular expressions matching the standard SDF property tag format:

```python
pattern = rf'>\s*<{field}>\s*\n([^\n]+)'
scores[field] = re.findall(pattern, content)
```

Extracted values are assembled into a pandas DataFrame, numeric columns cast with `pd.to_numeric(errors="coerce")`, and saved as CSV. The best pose is identified by minimum `minimizedAffinity` and converted from compressed SDF to PDB using Open Babel:

```bash
obabel kras_sotorasib_gnina_docked.sdf.gz -O kras_sotorasib_best_pose.pdb -f 1 -l 1
```

**Visualisations generated**

| Figure | Description | Output file |
|---|---|---|
| Affinity bar chart | `minimizedAffinity` for all 10 poses; y-axis inverted so stronger binding appears taller | `gnina_docking_scores.png` |
| Ligand-only 3D view | Best pose rendered in py3Dmol as magenta sticks, zoomed to the ligand | Inline Colab viewer |
| Receptor–ligand complex | Full KRAS receptor as spectrum cartoon + best pose Sotorasib as magenta sticks | `kras_sotorasib_receptor_ligand_visualization.png` |
| 3D scatter plot | Matplotlib 3D projection of receptor Cα backbone (line) and docked ligand atoms (scatter) | `kras_sotorasib_receptor_ligand_visualization.png` |

**Best pose final summary**

| Metric | Value |
|---|---|
| `minimizedAffinity` | −10.67 kcal/mol |
| `CNNscore` | 0.3321 (highest across all 10 poses) |
| `CNNaffinity` | 7.251 (highest across all 10 poses) |
| Intramolecular strain | 1.08 kcal/mol (acceptable) |

---

## 5. Clinical Relevance

Sotorasib (Lumakras) achieved a ~37% objective response rate in KRAS G12C-mutant NSCLC in the **CodeBreaK 100 Phase II trial** (*New England Journal of Medicine*, 2021). Its mechanism involves selective covalent modification of the mutant cysteine 12, stabilising KRAS in its inactive GDP-bound state and blocking downstream oncogenic signalling.

The predicted affinity of **−10.67 kcal/mol** and CNNscore of **0.3321** from this pipeline are consistent with Sotorasib's known potency (IC₅₀ ~0.18 nM in KRAS G12C biochemical assays), providing computational support for its experimentally validated binding mode.

---

## 6. Limitations

| Limitation | Detail |
|---|---|
| Rigid receptor docking | GNINA does not model backbone flexibility. KRAS switch I/II loops are highly dynamic — especially relevant to its activation mechanism. |
| Covalent inhibition not modelled | Sotorasib forms a permanent covalent bond with C12. Non-covalent docking cannot represent this; scores reflect binding site recognition only. |
| Wild-type sequence retrieved | `NM_004985.5` is the wild-type KRAS. The G12C mutation exists only in the 6OIM experimental structure — intentional by design. |
| Online BLAST | `NCBIWWW.qblast` depends on NCBI server load and may take 5–15 minutes. |
| AlphaFold single conformation | Predicts the most probable apo state only; does not capture ensemble dynamics or ligand-induced conformational changes. |

---

## 7. Key References

| Reference | Role in pipeline |
|---|---|
| McNutt AT et al. *J Cheminformatics* 2021;13:43 | GNINA 1.0 — docking engine used |
| Jumper J et al. *Nature* 2021;596:583–589 | AlphaFold — source of predicted structure |
| Skoulidis F et al. *NEJM* 2021;384:2371–2381 | CodeBreaK 100 — Sotorasib clinical trial |
| Hallin J et al. *Cancer Discov* 2020;10:54–71 | KRAS G12C biology and S-IIP mechanism |
| Cock PJA et al. *Bioinformatics* 2009;25:1422 | Biopython — sequence and structure handling |
| O'Boyle NM et al. *J Cheminformatics* 2011;3:33 | Open Babel — ligand format conversion |

---

*KRAS G12C Bioinformatics Pipeline — Methodology Report · May 2026*
