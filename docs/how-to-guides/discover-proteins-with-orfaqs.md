# Discover Proteins with ORFaqs
ORFaqs Protein Discovery processes genomic sequences and provides you with a
list of all possible proteins. The tool does models some biological aspects
of the converting genomic sequences into representative proteins
(transcription, translation, ect.). However, it makes no effort to do so in a
biologically plausible fashion.

Instead, the tool takes an unrestricted approach to protein discovery and
relies on two basic parameters:
- descriptions of the start and stop codons of interest, and
- the reading frames in which translation occurs.

With this approach, you are not constrained by what is *known*, but rather, you
are encouraged to ask: *What is possible?*

## Using the CLI
**Prerequisites**

If you are using the CLI from source, make sure you complete the items below before proceeding with this tutorial.
- [Setting Up Your Dev Environment](./../quick-start-quides/setting-up-your-dev-environment.md)
----

The CLI allows you to:
- Input one or more genomic sequences in FASTA file format.
- Input a single genomic sequence as a contiguous, unbroken string.

### Using FASTA Files as Input
```
orfaqs_protein_discovery -i input_fasta.fasta
```
<details>
<summary><i>Example FASTA File Format</i></summary>

Create a FASTA file by copying and pasting the text below into an empty file.
```
>YAL068C PAU8 SGDID:S000002142, Chr I from 2169-1807, Genome Release 64-5-1, reverse complement, Verified ORF, "Protein of unknown function; member of the seripauperin multigene family encoded mainly in subtelomeric regions"
ATGGTCAAATTAACTTCAATCGCCGCTGGTGTCGCTGCCATCGCTGCTACTGCTTCTGCA
ACCACCACTCTAGCTCAATCTGACGAAAGAGTCAACTTGGTGGAATTGGGTGTCTACGTC
TCTGATATCAGAGCTCACTTAGCCCAATACTACATGTTCCAAGCCGCCCACCCAACTGAA
ACCTACCCAGTCGAAGTTGCTGAAGCCGTTTTCAACTACGGTGACTTCACCACCATGTTG
ACCGGTATTGCTCCAGACCAAGTGACCAGAATGATCACCGGTGTTCCATGGTACTCCAGC
AGATTAAAGCCAGCCATCTCCAGTGCTCTATCCAAGGACGGTATCTACACTATCGCAAAC
TAG
```
</details>

### Using a Sequence String as Input

```
orfaqs_protein_discovery -i ATGGTCAAATTAACTTCAATCGCCGCTGGTGTCGCTGCCATCGCTGCTACTGCTTCTGCA
```

### Specifying the Exported Results Format
Three options are supported for exporting the discovered proteins:

| Export Option | Export File Extension |
|-|-|
| `--export_as_csv` *(default option)* |`.csv` |
| `--export_as_json` | `.json` |
| `--export_as_excel` | `.xlsx` |

#### Export Results in Excel Format
```
orfaqs_protein_discovery -i input_fasta.fasta --export_as_excel
```
#### Export Results in JSON Format
```
orfaqs_protein_discovery -i input_fasta.fasta --export_as_json
```

### Changing the Output Directory
All outputs from the app are generated within the directory path
`<OUTPUT_DIRECTORY>/.orfaqs-apps/orfaqs-protein-discovery`. By default,
`<OUTPUT_DIRECTORY>`is the current working directory.

Use the `-o, --output_directory` options to change the output directory path.
```
orfaqs_protein_discovery -i input_fasta.fasta -o /~/my-output-directory
```

### Organizing Results Using a Job ID
You can further organize your results by using a *job id*. The job id is joined
to the end of the output directory path resulting in the final path structure:
`<OUTPUT_DIRECTORY>/.orfaqs-apps/orfaqs-protein-discovery/<JOB_ID>`

Use the `-j, --job_id` options to add a job id directory to the output path.
```
orfaqs_protein_discovery -i input_fasta.fasta -o /~/my-output-directory -j my-job-id
```

## Output Results
### Directory Structure
As previously discussed, all outputs are placed in the
directory path:
`<OUTPUT_DIRECTORY>/.orfaqs-apps/orfaqs-protein-discovery/<JOB_ID>`.

<details>
<summary><i>(<b>NOTE</b>: Jobs ids are optional</i>)</summary>
If <code>-j, --job_id</code> options are not specified, a <code>JOB_ID</code> sub-directory will not appear as part of the output path.
</details>
</br>

Within this folder:
- **Outputs generated from FASTA files** appear in the path:
`<OUTPUT_DIRECTORY>/.orfaqs-apps/orfaqs-protein-discovery/<JOB_ID>/<SEQUENCE_ID>`,
where `<SEQUENCE_ID>` is the sequence identifier found in the FASTA file for the
corresponding sequence.
- **Outputs from sequence strings** appear directly in the directory path
`<OUTPUT_DIRECTORY>/<JOB_ID>/.orfaqs-apps/orfaqs-protein-discovery`.

### Output Files
Outputs are presented in two ways:
1. By reading frame
1. As a complete set of results (data from all reading frames included)

#### Outputs by Reading Frame
Each reading frame specific output file name is formatted as:
`discovered-proteins-reading-frame-<READING_FRAME>.<EXPORT_EXTENSION>` where:
- `<READING_FRAME>`is the index 1, 2, or 3, corresponding to the reading
frame used when the protein was found, and
- `<EXPORT_EXTENSION>` corresponds to the extension used by the export option selected (`.csv`, `.json`, or `.xlsx`)

If no proteins were found within a particular reading frame, then no output
data are exported for that reading frame.

#### All Discovered Proteins Outputs
Results for all the proteins discoverd for a given sequence are formatted as:
`discovered-proteins.<EXPORT_EXTENSION>` where:
- `<EXPORT_EXTENSION>` corresponds to the extension used by the export option selected (`.csv`, `.json`, or `.xlsx`)

If no proteins were found for a given protein sequence, then no output data are
exported.

#### Output File Contents
Within each file, results are organized either as tables (in the case of `CSV`
and `Excel` exports), or as lists of `JSON` objects (in the case of `JSON`
exports). Each file contains the following information:

- `reading_frame`: `(int)` The reading frame used during translation.
- `rna_sequence_position`: `(int)` The place in the RNA sequence where
translation began.
- `protein`: `(str)` The resulting protein sequence written using standard
single-letter code names.
- `protein_length`: `(int)` The total number of amino acids making up the
protein sequence.

The keys used in either format (columns for exported tables and keys in JSON
key-value pairs) are identical.
