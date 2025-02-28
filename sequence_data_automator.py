import pandas as pd
import requests
from collections import Counter
from itertools import product

def fetch_mirna_sequence(mirna_name):
    """Fetch miRNA sequence from miRBase (or another database if API available)."""
    base_url = "C:\Users\brand\Downloads\miRNA.dat"  # Example URL (adjust if needed)
    response = requests.get(base_url + mirna_name)
    
    return response.text.strip() if response.status_code == 200 else None

def generate_possible_motifs(length):
    """Generate all possible motifs of a given length (AA, AU, ... GGGG)."""
    return ["".join(m) for m in product("AUCG", repeat=length)]

def calculate_symmetry(seq):
    """Calculate symmetry score based on matching bases in the first half vs. reversed second half."""
    half = len(seq) // 2
    first_half, second_half = seq[:half], seq[-half:][::-1]
    return sum(1 for i in range(half) if first_half[i] == second_half[i])

def check_motif_frequencies(sequence, motifs):
    """Return a frequency vector for motifs in the sequence."""
    motif_counts = {motif: 0 for motif in motifs}
    
    for i in range(len(sequence) - len(motifs[0]) + 1):
        motif = sequence[i:i+len(motifs[0])]
        if motif in motif_counts:
            motif_counts[motif] += 1
    
    return [motif_counts[m] for m in motifs]

def calculate_sequence_descriptors(sequence):
    """Calculate various sequence descriptors for a given miRNA sequence."""
    counts = Counter(sequence)
    A, U, C, G = counts.get('A', 0), counts.get('U', 0), counts.get('C', 0), counts.get('G', 0)
    total = A + U + C + G

    fA, fU, fC, fG = (A/total, U/total, C/total, G/total) if total > 0 else (0, 0, 0, 0)
    mean_mass = (135.1*A + 112.1*U + 111.1*C + 151.1*G) / total if total > 0 else 0
    hydrogen_bonds = 2*(A + U) + 3*(C + G)
    symmetry_score = calculate_symmetry(sequence)

    # Generate all possible motifs
    motifs_2bp = generate_possible_motifs(2)
    motifs_3bp = generate_possible_motifs(3)
    motifs_4bp = generate_possible_motifs(4)

    # Check motif frequencies
    motif_2bp_values = check_motif_frequencies(sequence, motifs_2bp)
    motif_3bp_values = check_motif_frequencies(sequence, motifs_3bp)
    motif_4bp_values = check_motif_frequencies(sequence, motifs_4bp)

    return [A, U, C, G, fA, fU, fC, fG, mean_mass, hydrogen_bonds, symmetry_score] + motif_2bp_values + motif_3bp_values + motif_4bp_values

# Load miRNA list
mirna_list = [
    "hsa-miR-138-5p", "hsa-miR-338-3p", "hsa-miR-431-5p", "hsa-miR-30c-2-3p",
    "hsa-miR-10394-5p", "hsa-miR-1208", "hsa-miR-1256", "hsa-miR-1271-5p",
    "hsa-miR-1286", "hsa-miR-142-3p", "hsa-miR-183-5p", "hsa-miR-20a-3p",
    "hsa-miR-2276-3p", "hsa-miR-2278", "hsa-miR-27b-5p", "hsa-miR-299-3p",
    "hsa-miR-3158-5p", "hsa-miR-342-3p", "hsa-miR-3606-3p", "hsa-miR-3622a-3p",
    "hsa-miR-365a-3p", "hsa-miR-423-5p", "hsa-miR-4266", "hsa-miR-4272",
    "hsa-miR-449c-3p", "hsa-miR-452-3p", "hsa-miR-4660", "hsa-miR-4712-5p",
    "hsa-miR-4714-5p", "hsa-miR-4740-3p", "hsa-miR-4762-5p", "hsa-miR-4769-3p",
    "hsa-miR-496", "hsa-miR-5009-5p", "hsa-miR-526b-5p", "hsa-miR-548a-3p",
    "hsa-miR-574-3p", "hsa-miR-590-3p", "hsa-miR-6082", "hsa-miR-6128",
    "hsa-miR-627-5p", "hsa-miR-631", "hsa-miR-647", "hsa-miR-6719-3p",
    "hsa-miR-6727-5p", "hsa-miR-6757-3p", "hsa-miR-6779-3p", "hsa-miR-6792-3p",
    "hsa-miR-6796-5p", "hsa-miR-6811-3p", "hsa-miR-6848-3p", "hsa-miR-6873-3p",
    "hsa-miR-6874-3p", "hsa-miR-6882-5p", "hsa-miR-6889-3p", "hsa-miR-7162-5p"
]


# Generate column headers for motif presence
motif_2bp_cols = generate_possible_motifs(2)
motif_3bp_cols = generate_possible_motifs(3)
motif_4bp_cols = generate_possible_motifs(4)

data = []
for mirna in mirna_list:
    sequence = fetch_mirna_sequence(mirna)
    if sequence:
        descriptors = calculate_sequence_descriptors(sequence)
        data.append([mirna, sequence] + descriptors)
    else:
        # Append missing values (including symmetry score and motif presence)
        data.append([mirna, "Not Found"] + [""] * (11 + len(motif_2bp_cols) + len(motif_3bp_cols) + len(motif_4bp_cols)))

# Save to Excel
columns = (["miRNA Name", "Sequence", "A Count", "U Count", "C Count", "G Count", 
           "fA", "fU", "fC", "fG", "Mean Mass", "Hydrogen Bonds", "Symmetry Score"] 
           + motif_2bp_cols + motif_3bp_cols + motif_4bp_cols)

df = pd.DataFrame(data, columns=columns)
df.to_excel("miRNA_Descriptors3.xlsx", index=False)

print("Data saved to miRNA_Descriptors3.xlsx")
