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

kernel void find_codons(
    device packed_uchar3* rna_sequence [[buffer(0)]],
    constant uint &number_codons [[buffer(1)]],
    constant uint &data_stride [[buffer(2)]],
    device packed_uchar3* reference_codons [[buffer(3)]],
    constant uint &number_reference_codons [[buffer(4)]],
    device bool* found_codon_indices [[buffer(5)]],
    uint g_id [[thread_position_in_grid]]) {
    // Process all codons within the RNA sequence.
    for (uint g_i = g_id; g_i < number_codons; g_i += data_stride) {
        uchar3 current_codon = rna_sequence[g_i];
        for (int i = 0; i < number_reference_codons; ++i) {
            if (all(current_codon == reference_codons[i])) {
                found_codon_indices[g_i] = true;
                break;
            }
        }
    }
}

kernel void all_orf_lengths(
    device uint* start_codon_indices [[buffer(0)]],
    constant uint &number_start_codons [[buffer(1)]],
    device uint* stop_codon_indices [[buffer(2)]],
    constant uint &number_stop_codons [[buffer(3)]],
    constant uint &data_stride [[buffer(4)]],
    device uint* orf_lengths [[buffer(5)]],
    uint g_id [[thread_position_in_grid]]) {
    // Find the start and stop codon pair (if it exists)
    // defining the reading frame for the protein.
    for (int i = g_id; i < number_start_codons; i += data_stride) {
        int start_codon_index = start_codon_indices[i];
        int stop_codon_index = -1;
        bool start_codon_begins_orf = false;
        for (int j = 0; j < number_stop_codons; ++j) {
            stop_codon_index = stop_codon_indices[j];
            if (start_codon_index < stop_codon_index) {
                start_codon_begins_orf = true;
                break;
            }
        }

        // If the start codon marks the start of an ORF,
        // record the length of the ORF.
        if (start_codon_begins_orf) {
            // Do not include the stop codon when computing the length.
            orf_lengths[i] = stop_codon_index - start_codon_index;
        }
    }
}

kernel void rna_to_amino_acid_sequence(
    device packed_uchar3* rna_sequence [[buffer(0)]],
    constant uint &number_codons [[buffer(1)]],
    constant uint &data_stride [[buffer(2)]],
    device uchar* amino_acid_sequence [[buffer(3)]],
    uint g_id [[thread_position_in_grid]]) {
    for (int i = g_id; i < number_codons; i += data_stride) {
        amino_acid_sequence[i] = codon_to_amino_acid(rna_sequence[i]);
    }
}

kernel void translate_all_orfs(
    device uchar* amino_acid_sequence [[buffer(0)]],
    device uint* start_codon_indices [[buffer(1)]],
    device uint* orf_lengths [[buffer(2)]],
    constant uint &number_start_codons [[buffer(3)]],
    constant uint &data_stride [[buffer(4)]],
    device uchar* all_orf_proteins [[buffer(5)]],
    uint g_id [[thread_position_in_grid]]) {

    ////////////////////////////////////////////////////////////////////////////
    // Translate the ORF.
    for (int i = g_id; i < number_start_codons; i += data_stride) {
        // Determine the offset position within the protein buffer where the
        // resulting protein should be stored.
        uint protein_buffer_offset = 0;
        for (int j = 0; j < i; ++j) {
            protein_buffer_offset += orf_lengths[j];
        }

        // Find the start codon defining the ORF for the protein and copy the
        // pre-computed translations into the final output buffer.
        uint start_codon_index = start_codon_indices[i];
        device uchar* protein = &all_orf_proteins[protein_buffer_offset];
        uint protein_length = orf_lengths[i];
        for (int k = 0; k < protein_length; ++k) {
            protein[k] = amino_acid_sequence[start_codon_index + k];
        }
    }
}
