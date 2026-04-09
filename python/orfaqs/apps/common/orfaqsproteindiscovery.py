"""
ORFaqs Protein Discovery common app classes, resources, and utility functions.
"""

import enum
import logging
import pathlib

from tqdm import tqdm

from orfaqs.apps.common.orfaqsrecords import (
    ORFaqsDiscoveredProteinRecord,
    ORFaqsRecordUtils,
)

from orfaqs.lib.core.codons import Codon
from orfaqs.lib.core.enzymes import RNAPolymerase
from orfaqs.lib.utils.fastautils import FASTAUtils
from orfaqs.lib.core.nucleotides import (
    DNASequence,
    StrandType,
    GenomicSequence,
    NucleotideUtils,
    RNASequence,
)
from orfaqs.lib.core.proteins import Protein
from orfaqs.lib.core.ribosomes import (
    RibosomeUtils,
    RNAReadingFrame,
)

from orfaqs.lib.utils.directoryutils import DirectoryUtils
from orfaqs.lib.utils.pandasutils import (
    DataFrameExportFormat,
    DataFrameExportFormatOptions,
    PandasUtils,
)
from orfaqs.lib.utils.perfutils import PerfProfiler


_logger = logging.getLogger(__name__)
_perf_profiler = PerfProfiler(__name__)


class _ProfilingFunctionName(enum.Enum):
    PROTEIN_DISCOVERY_GENOMIC_SEQUENCE = enum.auto()
    TRANSLATE_RNA_GROUP = enum.auto()


_ExportFormatOptions = DataFrameExportFormatOptions
_AVAILABLE_EXPORT_FORMATS: list[str] = [
    format for format in _ExportFormatOptions.__args__
]


class ORFaqsProteinDiscoveryApi:
    """ORFaqsProteinDiscoveryApi"""

    DATAFRAME_INDEX_KEY = 'index'

    @staticmethod
    def default_export_format() -> str:
        return DataFrameExportFormat.CSV

    @staticmethod
    def default_output_directory() -> str:
        return './'

    @staticmethod
    def available_export_formats() -> list[str]:
        return _AVAILABLE_EXPORT_FORMATS

    @staticmethod
    def exported_dataframe_keys() -> list[str]:
        return [
            ORFaqsProteinDiscoveryApi.DATAFRAME_INDEX_KEY
        ] + ORFaqsDiscoveredProteinRecord.keys()

    @staticmethod
    def unreferenced_uid() -> str:
        "unknown_reference"

    @staticmethod
    def _exported_discovered_proteins_file_name(
        export_format: _ExportFormatOptions,
        uid: str = None,
    ) -> str:
        file_name = f'discovered-proteins.{export_format}'
        if uid is not None:
            file_name = f'{uid}-{file_name}'
        return file_name

    @staticmethod
    def _export_discovered_proteins(
        file_path: str | pathlib.Path,
        discovered_proteins: list[ORFaqsDiscoveredProteinRecord],
        export_format: _ExportFormatOptions,
    ):
        records_dataframe = ORFaqsRecordUtils.orfaqs_records_to_dataframe(
            discovered_proteins
        )
        file_path = DirectoryUtils.make_path_object(file_path)
        PandasUtils.export_dataframe(
            file_path=file_path,
            dataframe=records_dataframe,
            export_format=export_format,
            include_index=True,
        )

    @staticmethod
    def _discover_proteins(
        rna_sequence: RNASequence,
        reading_frames: list[RNAReadingFrame],
        start_codons: list[Codon],
        stop_codons: list[Codon],
        use_gpu: bool,
    ) -> dict[RNAReadingFrame, dict[int, Protein]]:
        #######################################################################
        # Initialize the protein map and related variables which store the
        # the results of the discovered proteins.
        reading_frame_proteins_map: dict[
            RNAReadingFrame, dict[int, Protein]
        ] = {}
        for reading_frame in reading_frames:
            (start_index, stop_index) = (
                RibosomeUtils.sequence_start_stop_indices(
                    rna_sequence, reading_frame
                )
            )
            ###################################################################
            # Get all proteins for the current reading frame.
            rna_sequence_frame = rna_sequence[start_index:stop_index]
            reading_frame_proteins_map[reading_frame] = (
                RibosomeUtils.translate_all_orfs(
                    rna_sequence=rna_sequence_frame,
                    start_codons=start_codons,
                    stop_codons=stop_codons,
                    use_gpu=use_gpu,
                )
            )
        return reading_frame_proteins_map

    @staticmethod
    def discover_proteins(
        genomic_sequence: str | GenomicSequence,
        uid: str = None,
        strand_type: StrandType = None,
        frames: list[RNAReadingFrame] = None,
        start_codons: list[Codon] = None,
        stop_codons: list[Codon] = None,
        include_reverse_complement: bool = True,
        output_directory: str | pathlib.Path = None,
        export_format: _ExportFormatOptions = None,
        use_gpu: bool = True,
    ) -> list[ORFaqsDiscoveredProteinRecord]:
        _perf_profiler.start_perf_timer(
            _ProfilingFunctionName.PROTEIN_DISCOVERY_GENOMIC_SEQUENCE
        )
        output_directory = DirectoryUtils.make_path_object(output_directory)
        DirectoryUtils.mkdir_path(output_directory)

        rna_sequence: RNASequence = None
        if isinstance(genomic_sequence, DNASequence):
            rna_sequence = RNAPolymerase.transcribe(genomic_sequence)
        elif isinstance(genomic_sequence, RNASequence):
            rna_sequence = genomic_sequence
        else:
            message = (
                '[ERROR] Expected type GenomicSequence.\n'
                '(debug) ->\n'
                f'\ttype(genomic_sequence): {type(genomic_sequence)}'
            )
            raise ValueError(message)

        if uid is None:
            uid = ORFaqsProteinDiscoveryApi.unreferenced_uid()
        if not isinstance(rna_sequence, RNASequence):
            rna_sequence = RNASequence(rna_sequence, strand_type=strand_type)

        reading_frames: list[RNAReadingFrame] = []
        if frames is None:
            reading_frames = RibosomeUtils.available_reading_frames()
        elif isinstance(frames, (list, set)):
            reading_frames = list(set(frames))
        elif isinstance(frames, RNAReadingFrame):
            reading_frames = [frames]

        strand_type_protein_map: dict[
            StrandType, dict[RNAReadingFrame, dict[int, Protein]]
        ] = {}
        #######################################################################
        # 1. Discover proteins for each reading frame in the give strand_type
        # direction.
        strand_type_protein_map[rna_sequence.strand_type] = (
            ORFaqsProteinDiscoveryApi._discover_proteins(
                rna_sequence=rna_sequence,
                reading_frames=reading_frames,
                start_codons=start_codons,
                stop_codons=stop_codons,
                use_gpu=use_gpu,
            )
        )

        # Process the reverse complement of the sequence if it is requested.
        if include_reverse_complement:
            rna_sequence = NucleotideUtils.reverse_complement(rna_sequence)
            strand_type_protein_map[rna_sequence.strand_type] = (
                ORFaqsProteinDiscoveryApi._discover_proteins(
                    rna_sequence=rna_sequence,
                    reading_frames=reading_frames,
                    start_codons=start_codons,
                    stop_codons=stop_codons,
                    use_gpu=use_gpu,
                )
            )

        _perf_profiler.stop_perf_timer(
            _ProfilingFunctionName.PROTEIN_DISCOVERY_GENOMIC_SEQUENCE
        )
        discovered_proteins: list[ORFaqsDiscoveredProteinRecord] = []
        for (
            strand_type,
            reading_frames_proteins_map,
        ) in strand_type_protein_map.items():
            for (
                reading_frame,
                proteins_map,
            ) in reading_frames_proteins_map.items():
                for base_index, protein in proteins_map.items():
                    # Convert the protein's base location to a "1's" based index.
                    discovered_proteins.append(
                        ORFaqsDiscoveredProteinRecord(
                            uid=uid,
                            strand_type=strand_type,
                            reading_frame=reading_frame,
                            rna_sequence_position=(base_index + 1),
                            protein=protein,
                        )
                    )

        if export_format is not None:
            file_name = ORFaqsProteinDiscoveryApi._exported_discovered_proteins_file_name(
                export_format=export_format,
                uid=uid,
            )
            file_path = output_directory / file_name
            ORFaqsProteinDiscoveryApi._export_discovered_proteins(
                file_path=file_path,
                discovered_proteins=discovered_proteins,
                export_format=export_format,
            )
        return discovered_proteins

    @staticmethod
    def process_genomic_sequence(
        genomic_sequence: str | GenomicSequence,
        strand_type: StrandType = None,
        uid: str = None,
        include_reverse_complement: bool = True,
        export_format: _ExportFormatOptions = None,
        output_directory: str | pathlib.Path = None,
        use_gpu: bool = True,
    ) -> list[ORFaqsDiscoveredProteinRecord]:
        if export_format is None:
            export_format = ORFaqsProteinDiscoveryApi.default_export_format()

        if output_directory is None:
            output_directory = (
                ORFaqsProteinDiscoveryApi.default_output_directory()
            )

        if isinstance(genomic_sequence, str):
            genomic_sequence = NucleotideUtils.create_sequence(
                genomic_sequence,
                strand_type=strand_type,
            )

        return ORFaqsProteinDiscoveryApi.discover_proteins(
            genomic_sequence,
            uid=uid,
            include_reverse_complement=include_reverse_complement,
            output_directory=output_directory,
            export_format=export_format,
            use_gpu=use_gpu,
        )

    @staticmethod
    def process_fasta_file(
        fasta_file_path: str | pathlib.Path,
        strand_type: StrandType = None,
        include_reverse_complement: bool = True,
        export_format: _ExportFormatOptions = None,
        output_directory: str | pathlib.Path = None,
        display_progress: bool = True,
        use_gpu: bool = True,
    ) -> dict[str, list[ORFaqsDiscoveredProteinRecord]]:
        if export_format is None:
            export_format = ORFaqsProteinDiscoveryApi.default_export_format()

        if output_directory is None:
            output_directory = (
                ORFaqsProteinDiscoveryApi.default_output_directory()
            )

        fasta_sequences = FASTAUtils.parse_file(fasta_file_path)
        if output_directory is not None:
            output_directory = DirectoryUtils.make_path_object(
                output_directory
            )
        discovered_proteins_maps: dict[
            str, list[ORFaqsDiscoveredProteinRecord]
        ] = {}
        total_protein_count = 0
        process_list = fasta_sequences
        if display_progress:
            tqdm_description = 'Processing genomic sequences...'
            process_list = tqdm(
                process_list, desc=tqdm_description, total=len(process_list)
            )

        for fasta_sequence in process_list:
            print('-------------------------------------------------------')
            print(f'Processing Sequence: {fasta_sequence.name}')
            current_sequence_output_directory = output_directory.joinpath(
                f'{fasta_sequence.uid}'
            )
            discovered_proteins = (
                ORFaqsProteinDiscoveryApi.process_genomic_sequence(
                    fasta_sequence.sequence,
                    strand_type=strand_type,
                    uid=fasta_sequence.uid,
                    include_reverse_complement=include_reverse_complement,
                    output_directory=current_sequence_output_directory,
                    export_format=export_format,
                    use_gpu=use_gpu,
                )
            )

            discovered_proteins_maps[fasta_sequence.uid] = discovered_proteins
            total_protein_count += len(discovered_proteins)
            print('\n')

        print(f'Total number of proteins discovered: {total_protein_count}')
        return (discovered_proteins_maps, total_protein_count)
