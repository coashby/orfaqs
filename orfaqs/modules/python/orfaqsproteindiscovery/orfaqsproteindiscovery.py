"""
ORFaqs Protein Discovery common app classes, resources, and utility functions.
"""

import enum
import logging
import os
import pathlib

from tqdm import tqdm

from orfaqs.modules.python.orfaqsrecords.orfaqsrecords import (
    ORFaqsDiscoveredProteinRecord,
    ORFaqsRecordUtils,
)

from orfaqs.lib.python.core.codons import Codon
from orfaqs.lib.python.core.enzymes import RNAPolymerase
from orfaqs.lib.python.utils.fastautils import FASTAUtils
from orfaqs.lib.python.core.nucleotides import (
    DNASequence,
    StrandType,
    GenomicSequence,
    NucleotideUtils,
    RNASequence,
)
from orfaqs.lib.python.core.proteins import Protein
from orfaqs.lib.python.core.ribosomes import (
    RibosomeUtils,
    RNAReadingFrame,
)

from orfaqs.lib.python.utils.directoryutils import DirectoryUtils
from orfaqs.lib.python.utils.pandasutils import (
    DataFrameExportFormat,
    DataFrameExportFormatOptions,
    PandasUtils,
)
from orfaqs.lib.python.utils.perfutils import PerfProfiler


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
        discovered_proteins: list[ORFaqsDiscoveredProteinRecord],
        output_directory: str | os.PathLike = None,
        uid: str = None,
        export_format: _ExportFormatOptions = None,
    ):
        records_dataframe = ORFaqsRecordUtils.orfaqs_records_to_dataframe(
            discovered_proteins
        )

        if output_directory is None:
            output_directory = (
                ORFaqsProteinDiscoveryApi.default_output_directory()
            )

        if export_format is None:
            export_format = ORFaqsProteinDiscoveryApi.default_export_format()

        output_directory = DirectoryUtils.make_path_object(output_directory)
        DirectoryUtils.mkdir_path(output_directory)

        file_name = (
            ORFaqsProteinDiscoveryApi._exported_discovered_proteins_file_name(
                export_format=export_format,
                uid=uid,
            )
        )
        file_path = DirectoryUtils.make_path_object(
            output_directory / file_name
        )
        PandasUtils.export_dataframe(
            file_path=file_path,
            dataframe=records_dataframe,
            export_format=export_format,
            include_index=True,
        )

        return file_path

    @staticmethod
    def _translate_all_orf(
        rna_sequence: RNASequence,
        reading_frames: list[RNAReadingFrame],
        start_codons: list[Codon],
        stop_codons: list[Codon],
        use_gpu: bool,
        display_progress: bool = False,
    ) -> dict[RNAReadingFrame, dict[int, Protein]]:
        #######################################################################
        # Initialize the protein map and related variables which store the
        # the results of the discovered proteins.
        reading_frame_proteins_map: dict[
            RNAReadingFrame, dict[int, Protein]
        ] = {}
        for reading_frame in reading_frames:
            if display_progress:
                print(f'Reading Frame: {reading_frame}')
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
                    display_progress=display_progress,
                )
            )
        return reading_frame_proteins_map

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
        return 'unknown_reference'

    @staticmethod
    def find_discovered_protein_files(
        input_path: str | os.PathLike,
    ) -> list[pathlib.Path]:
        #######################################################################
        # Gather all discovery files.
        input_file_paths: list[pathlib.Path] = []
        if DirectoryUtils.is_file(input_path):
            input_file_paths.append(input_path)

        elif DirectoryUtils.is_directory(input_path):
            # Grab all files from the directory. In discoveries involving
            # multiple genes, discovery files will be organized in
            # subdirectories.
            input_file_paths = DirectoryUtils.glob_files(
                input_path, recursive=True
            )

        discovered_proteins_files: list[pathlib.Path] = []
        for file_path in input_file_paths:
            # Validate the file paths...
            proteins_dataframe = PandasUtils.read_file_as_dataframe(
                file_path,
                raise_error=False,
            )
            if proteins_dataframe is None:
                continue

            for record_key in ORFaqsDiscoveredProteinRecord.keys():
                if record_key not in proteins_dataframe.columns:
                    continue
            discovered_proteins_files.append(file_path)

        return discovered_proteins_files

    @staticmethod
    def discover_proteins_from_sequence(
        genomic_sequence: str | GenomicSequence,
        uid: str = None,
        strand_type: StrandType = None,
        frames: list[RNAReadingFrame] = None,
        start_codons: list[Codon] = None,
        stop_codons: list[Codon] = None,
        include_reverse_complement: bool = True,
        export_results: bool = False,
        output_directory: str | os.PathLike = None,
        export_format: _ExportFormatOptions = None,
        use_gpu: bool = True,
        display_progress: bool = False,
    ) -> tuple[(list[ORFaqsDiscoveredProteinRecord] | pathlib.Path), int]:
        _perf_profiler.start_perf_timer(
            _ProfilingFunctionName.PROTEIN_DISCOVERY_GENOMIC_SEQUENCE
        )
        if isinstance(genomic_sequence, str):
            genomic_sequence = NucleotideUtils.make_sequence_object(
                genomic_sequence,
                strand_type=strand_type,
            )

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
            ORFaqsProteinDiscoveryApi._translate_all_orf(
                rna_sequence=rna_sequence,
                reading_frames=reading_frames,
                start_codons=start_codons,
                stop_codons=stop_codons,
                use_gpu=use_gpu,
                display_progress=display_progress,
            )
        )

        # Process the reverse complement of the sequence if it is requested.
        if include_reverse_complement:
            rna_sequence = NucleotideUtils.reverse_complement(rna_sequence)
            strand_type_protein_map[rna_sequence.strand_type] = (
                ORFaqsProteinDiscoveryApi._translate_all_orf(
                    rna_sequence=rna_sequence,
                    reading_frames=reading_frames,
                    start_codons=start_codons,
                    stop_codons=stop_codons,
                    use_gpu=use_gpu,
                    display_progress=display_progress,
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
                if proteins_map is None:
                    continue
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
        number_proteins = len(discovered_proteins)
        if export_results or export_format is not None:
            file_path = ORFaqsProteinDiscoveryApi._export_discovered_proteins(
                discovered_proteins=discovered_proteins,
                output_directory=output_directory,
                uid=uid,
                export_format=export_format,
            )
            return (file_path, number_proteins)

        return (discovered_proteins, number_proteins)

    @staticmethod
    def discover_proteins_from_fasta_file(
        fasta_file_path: str | os.PathLike,
        strand_type: StrandType = None,
        include_reverse_complement: bool = True,
        export_results: bool = False,
        output_directory: str | os.PathLike = None,
        export_format: _ExportFormatOptions = None,
        display_progress: bool = False,
        use_gpu: bool = True,
    ) -> tuple[dict[str, list[ORFaqsDiscoveredProteinRecord]], int]:
        fasta_sequences = FASTAUtils.parse_file(fasta_file_path)
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

        # Create any required output paths.
        if export_format is not None:
            export_results = True

        if export_results and output_directory is None:
            output_directory = (
                ORFaqsProteinDiscoveryApi.default_output_directory()
            )
        if export_results:
            output_directory = DirectoryUtils.make_path_object(
                output_directory
            )
        current_sequence_output_directory = None
        for fasta_sequence in process_list:
            if display_progress:
                print()
                print(
                    '-------------------------------------------------------'
                )
                print(f'Processing Sequence: {fasta_sequence.name}')
            if export_results:
                current_sequence_output_directory = output_directory.joinpath(
                    f'{fasta_sequence.uid}'
                )
            (discovered_proteins, number_proteins) = (
                ORFaqsProteinDiscoveryApi.discover_proteins_from_sequence(
                    fasta_sequence.sequence,
                    strand_type=strand_type,
                    uid=fasta_sequence.uid,
                    include_reverse_complement=include_reverse_complement,
                    export_results=export_results,
                    output_directory=current_sequence_output_directory,
                    export_format=export_format,
                    use_gpu=use_gpu,
                    display_progress=display_progress,
                )
            )

            discovered_proteins_maps[fasta_sequence.uid] = discovered_proteins
            total_protein_count += number_proteins

        if display_progress:
            print(
                f'Total number of proteins discovered: {total_protein_count}'
            )
        return (discovered_proteins_maps, total_protein_count)

    @staticmethod
    def discover_proteins(
        input_sequence: str | os.PathLike,
        uid: str = None,
        include_reverse_complement: bool = True,
        output_directory: str | os.PathLike = None,
        job_id: str = None,
        export_results: bool = False,
        export_format: str = None,
        enable_gpu: bool = False,
        display_progress: bool = False,
    ) -> tuple[any, int]:
        #######################################################################
        # Create the local output directory path
        if output_directory is None:
            output_directory = './'

        output_directory = DirectoryUtils.make_path_object(output_directory)
        output_directory = output_directory.joinpath(
            ORFaqsProteinDiscoveryApi.default_output_directory()
        )
        output_directory = DirectoryUtils.make_path_object(output_directory)
        if isinstance(job_id, str):
            output_directory = output_directory.joinpath(job_id)

        DirectoryUtils.mkdir_path(output_directory)

        if DirectoryUtils.is_file(input_sequence):
            # Process as a FASTA file
            return ORFaqsProteinDiscoveryApi.discover_proteins_from_fasta_file(
                fasta_file_path=input_sequence,
                include_reverse_complement=include_reverse_complement,
                export_results=export_results,
                export_format=export_format,
                output_directory=output_directory,
                use_gpu=enable_gpu,
                display_progress=display_progress,
            )
        else:
            # Process as an sequence string.
            return ORFaqsProteinDiscoveryApi.discover_proteins_from_sequence(
                genomic_sequence=input_sequence,
                uid=uid,
                include_reverse_complement=include_reverse_complement,
                export_results=export_results,
                output_directory=output_directory,
                export_format=export_format,
                use_gpu=enable_gpu,
                display_progress=display_progress,
            )
