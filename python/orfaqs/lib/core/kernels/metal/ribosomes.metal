#include <metal_stdlib>
using namespace metal;

/*
Find all start codons in the given sequence and return the indices in the output
buffer.

:param rna_sequence packed_uchar3*:
    An RNA sequence stored as a packed buffer of uchar3 vectors.
:param start_codons packed_uchar3*:
    <>
:param number_start_codons const int:
    <>
:param start_codon_indices bool*:
    <>
*/
kernel void find_start_codons(
    device packed_uchar3* rna_sequence [[buffer(0)]],
    device packed_uchar3* start_codons [[buffer(1)]],
    device int* number_start_codons [[buffer(2)]],
    device bool* start_codon_indices [[buffer(3)]],
    uint g_id [[thread_position_in_grid]]) {
    uchar3 current_codon = rna_sequence[g_id];
    for (int i = 0; i < number_start_codons[0]; ++i) {
        if (all(current_codon == start_codons[i])) {
            start_codon_indices[g_id] = true;
            break;
        }
    }
}
