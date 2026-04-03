#include <metal_stdlib>
using namespace metal;


uchar codon_to_amino_acid(uchar3 codon) {
    uchar amino_acid;
    if (all(codon == uchar3('G', 'C', 'U')) ||
        all(codon == uchar3('G', 'C', 'C')) ||
        all(codon == uchar3('G', 'C', 'A')) ||
        all(codon == uchar3('G', 'C', 'G'))) {
        amino_acid = uchar('A');
    } else if (all(codon == uchar3('C', 'G', 'U')) ||
               all(codon == uchar3('C', 'G', 'C')) ||
               all(codon == uchar3('C', 'G', 'A')) ||
               all(codon == uchar3('C', 'G', 'G')) ||
               all(codon == uchar3('A', 'G', 'A')) ||
               all(codon == uchar3('A', 'G', 'G'))) {
        amino_acid = uchar('R');
    } else if (all(codon == uchar3('A', 'A', 'U')) ||
               all(codon == uchar3('A', 'A', 'C'))) {
            amino_acid = uchar('N');
    } else if (all(codon == uchar3('G', 'A', 'U')) ||
               all(codon == uchar3('G', 'A', 'C'))) {
            amino_acid = uchar('D');
    } else if (all(codon == uchar3('U', 'G', 'U')) ||
               all(codon == uchar3('U', 'G', 'C'))) {
            amino_acid = uchar('C');
    } else if (all(codon == uchar3('C', 'A', 'A')) ||
               all(codon == uchar3('C', 'A', 'G'))) {
            amino_acid = uchar('Q');
    } else if (all(codon == uchar3('G', 'A', 'A')) ||
               all(codon == uchar3('G', 'A', 'G'))) {
            amino_acid = uchar('E');
    } else if (all(codon == uchar3('G', 'G', 'U')) ||
               all(codon == uchar3('G', 'G', 'C')) ||
               all(codon == uchar3('G', 'G', 'A')) ||
               all(codon == uchar3('G', 'G', 'G'))) {
            amino_acid = uchar('G');
    } else if (all(codon == uchar3('C', 'A', 'U')) ||
               all(codon == uchar3('C', 'A', 'C'))) {
            amino_acid = uchar('H');
    } else if (all(codon == uchar3('A', 'U', 'U')) ||
               all(codon == uchar3('A', 'U', 'C')) ||
               all(codon == uchar3('A', 'U', 'A'))) {
            amino_acid = uchar('I');
    } else if (all(codon == uchar3('U', 'U', 'A')) ||
               all(codon == uchar3('U', 'U', 'G')) ||
               all(codon == uchar3('C', 'U', 'U')) ||
               all(codon == uchar3('C', 'U', 'C')) ||
               all(codon == uchar3('C', 'U', 'A')) ||
               all(codon == uchar3('C', 'U', 'G'))) {
            amino_acid = uchar('L');
    } else if (all(codon == uchar3('A', 'A', 'A')) ||
               all(codon == uchar3('A', 'A', 'G'))) {
            amino_acid = uchar('K');
    } else if (all(codon == uchar3('A', 'U', 'G'))) {
            amino_acid = uchar('M');
    } else if (all(codon == uchar3('U', 'U', 'U')) ||
               all(codon == uchar3('U', 'U', 'C'))) {
            amino_acid = uchar('F');
    } else if (all(codon == uchar3('C', 'C', 'U')) ||
               all(codon == uchar3('C', 'C', 'C')) ||
               all(codon == uchar3('C', 'C', 'A')) ||
               all(codon == uchar3('C', 'C', 'G'))) {
            amino_acid = uchar('P');
    } else if (all(codon == uchar3('U', 'C', 'U')) ||
               all(codon == uchar3('U', 'C', 'C')) ||
               all(codon == uchar3('U', 'C', 'A')) ||
               all(codon == uchar3('U', 'C', 'G')) ||
               all(codon == uchar3('A', 'G', 'U')) ||
               all(codon == uchar3('A', 'G', 'C'))) {
            amino_acid = uchar('S');
    } else if (all(codon == uchar3('A', 'C', 'U')) ||
               all(codon == uchar3('A', 'C', 'C')) ||
               all(codon == uchar3('A', 'C', 'A')) ||
               all(codon == uchar3('A', 'C', 'G'))) {
            amino_acid = uchar('T');
    } else if (all(codon == uchar3('U', 'G', 'G'))) {
            amino_acid = uchar('W');
    } else if (all(codon == uchar3('U', 'A', 'U')) ||
               all(codon == uchar3('U', 'A', 'C'))) {
            amino_acid = uchar('Y');
    } else if (all(codon == uchar3('G', 'U', 'U')) ||
               all(codon == uchar3('G', 'U', 'C')) ||
               all(codon == uchar3('G', 'U', 'A')) ||
               all(codon == uchar3('G', 'U', 'G'))) {
            amino_acid = uchar('V');
    } else if (all(codon == uchar3('U', 'G', 'A'))) {
            amino_acid = uchar('U');
    } else if (all(codon == uchar3('U', 'A', 'G'))) {
            amino_acid = uchar('O');
    }
    return amino_acid;
}

kernel void all_orf_lengths(
    device uint* start_codon_indices [[buffer(0)]],
    constant uint &number_start_codons [[buffer(1)]],
    device uint* stop_codon_indices [[buffer(2)]],
    constant uint &number_stop_codons [[buffer(3)]],
    device uint* orf_lengths [[buffer(4)]],
    uint g_id [[thread_position_in_grid]]) {
    // Find the start and stop codon pair (if it exists)
    // defining the reading frame for the protein.
    int start_codon_index = start_codon_indices[g_id];
    int stop_codon_index = -1;
    bool start_codon_begins_orf = false;
    for (int i = 0; i < number_stop_codons; ++i) {
        stop_codon_index = stop_codon_indices[i];
        if (start_codon_index < stop_codon_index) {
            start_codon_begins_orf = true;
            break;
        }
    }

    // If the start codon marks the start of an ORF,
    // record the length of the ORF.
    if (start_codon_begins_orf) {
        orf_lengths[g_id] = stop_codon_index - start_codon_index;
    }

}
kernel void translate_all_orfs(
    device packed_uchar3* rna_sequence [[buffer(0)]],
    device uint* start_codon_indices [[buffer(1)]],
    constant uint &number_start_codons [[buffer(2)]],
    device uint* orf_lengths [[buffer(3)]],
    device uchar* all_orf_proteins [[buffer(4)]],
    uint g_id [[thread_position_in_grid]]) {

    if (g_id >= number_start_codons) {
        return;
    }
    // Find the start and stop codon pair (if it exists)
    // defining the reading frame for the protein.
    int start_codon_index = start_codon_indices[g_id];

    // Translate the ORF.
    uint protein_buffer_offset = 0;
    for (int i = 0; i < g_id; ++i) {
        protein_buffer_offset += orf_lengths[i];
    }

    device uchar* protein = &all_orf_proteins[protein_buffer_offset];
    uint protein_length = orf_lengths[g_id];
    for (int i = 0; i < protein_length; ++i) {
        uchar3 codon = rna_sequence[start_codon_index + i];
        protein[i] = codon_to_amino_acid(codon);
    }
}

/*
Find all start codons in the given sequence and return the indices in the output
buffer.

:param rna_sequence packed_uchar3*:
    An RNA sequence stored as a packed buffer of uchar3 vectors.
:param codons packed_uchar3*:
    <>
:param number_start_codons const int:
    <>
:param start_codon_indices bool*:
    <>
*/
kernel void find_codons(
    device packed_uchar3* rna_sequence [[buffer(0)]],
    device packed_uchar3* codons [[buffer(1)]],
    constant uint &number_start_codons [[buffer(2)]],
    device bool* start_codon_indices [[buffer(3)]],
    uint g_id [[thread_position_in_grid]]) {
    uchar3 current_codon = rna_sequence[g_id];
    for (int i = 0; i < number_start_codons; ++i) {
        if (all(current_codon == codons[i])) {
            start_codon_indices[g_id] = true;
            break;
        }
    }
}
