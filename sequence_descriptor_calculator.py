from collections import Counter

def calculate_sequence_descriptors(sequence):
    # count string length
    total_length = len(sequence)
    
    # Count nucleotides
    counts = Counter(sequence)
    A, U, C, G = counts.get('A', 0), counts.get('U', 0), counts.get('C', 0), counts.get('G', 0)
    total = A + U + C + G

    # Frequency calculations
    fA, fU, fC, fG = (A/total, U/total, C/total, G/total) if total > 0 else (0, 0, 0, 0)

    # Mean mass
    mean_mass = (135.1*A + 112.1*U + 111.1*C + 151.1*G) / total if total > 0 else 0

    # Hydrogen bonds
    hydrogen_bonds = 2*(A + U) + 3*(C + G)

    # Symmetry Score: Count matching bases in the first half vs. reversed second half
    def calculate_symmetry(seq):
        length = len(seq)
        half = length // 2
        first_half = seq[:half]
        second_half = seq[-half:][::-1]  # Reverse second half
        matches = sum(1 for i in range(half) if first_half[i] == second_half[i])
        return matches  # Return as an integer count

    symmetry_score = calculate_symmetry(sequence)

    # Motif checks
    def check_motifs(sequence, length):
        motifs = {}
        for i in range(len(sequence) - length + 1):
            motif = sequence[i:i+length]
            motifs[motif] = 1
        return motifs

    motifs_2bp = check_motifs(sequence, 2)
    motifs_3bp = check_motifs(sequence, 3)
    motifs_4bp = check_motifs(sequence, 4)

    # Print descriptors only (no extra text)
    print(total_length)
    print(A, U, C, G)
    print(fA, fU, fC, fG)
    print(mean_mass)
    print(hydrogen_bonds)
    print(symmetry_score)  # Now an integer
    print(" ".join(motifs_2bp.keys()))
    print(" ".join(motifs_3bp.keys()))
    print(" ".join(motifs_4bp.keys()))

# Example usage with a sample sequence
sample_sequence = "AGCUGGUGUUGUGAAUCAGGCCG"
calculate_sequence_descriptors(sample_sequence)
