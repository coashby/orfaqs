"""
ORFaqs Protein Discovery common app classes, resources, and utility functions.
"""

import enum
import logging
import math
import multiprocessing
import os
import pandas as pd
import pathlib

from datetime import datetime
from tqdm import tqdm


from orfaqs.lib.core.codons import Codon
from orfaqs.lib.core.enzymes import RNAPolymerase
from orfaqs.lib.utils.fastautils import FASTAUtils
from orfaqs.lib.core.nucleotides import (
    DNASequence,
    GenomicSequence,
    NucleotideUtils,
    RNASequence,
)
from orfaqs.lib.core.proteins import Protein
from orfaqs.lib.core.ribosomes import (
    Ribosome,
    RibosomeUtils,
    RNAReadingFrame,
)

from orfaqs.lib.utils.directoryutils import DirectoryUtils
from orfaqs.lib.utils.jsonutils import JsonUtils
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


class ORFaqsProteinDiscoveryRecordKeys:
    """ORFaqsProteinDiscoveryRecordKeys"""

    ACCESSION_NUMBER_KEY = 'accession_number'
    PROTEIN_KEY = 'protein'
    PROTEIN_LENGTH_KEY = 'protein_length'
    READING_FRAME_KEY = 'reading_frame'
    RNA_SEQUENCE_POSITION_KEY = 'rna_sequence_position'


class ORFaqsProteinDiscoveryRecord(ORFaqsProteinDiscoveryRecordKeys):
    """ORFaqsProteinDiscoveryRecord"""

    def __init__(
        self,
        accession_number: str,
        reading_frame: RNAReadingFrame,
        rna_sequence_position: int,
        protein: Protein,
    ):
        self._accession_number = accession_number
        self._reading_frame = reading_frame
        self._rna_sequence_position = rna_sequence_position
        self._protein = protein

    @property
    def record(self) -> dict[str, any]:
        record_map: dict[str, any] = {}
        for record_key in ORFaqsProteinDiscoveryRecord.keys():
            if ORFaqsProteinDiscoveryRecord.ACCESSION_NUMBER_KEY == record_key:
                record_map[record_key] = self._accession_number
            elif ORFaqsProteinDiscoveryRecord.READING_FRAME_KEY == record_key:
                record_map[record_key] = self._reading_frame.value
            elif (
                ORFaqsProteinDiscoveryRecord.RNA_SEQUENCE_POSITION_KEY
                == record_key
            ):
                record_map[record_key] = self._rna_sequence_position
            elif ORFaqsProteinDiscoveryRecord.PROTEIN_KEY == record_key:
                record_map[record_key] = self._protein
            elif ORFaqsProteinDiscoveryRecord.PROTEIN_LENGTH_KEY == record_key:
                record_map[record_key] = self._protein.number_amino_acids
            else:
                message = (
                    '[ERROR] Missing record key assignment.\n'
                    '(debug) ->\n'
                    f'\trecord_key: {record_key} (not assigned)'
                )
                _logger.error(message)
                raise RuntimeError(message)

        return record_map

    @property
    def condensed_record_json_str(self) -> str:
        record_json_str = JsonUtils.as_json_string(
            JsonUtils.make_writable(self.record), indent=None
        )
        record_json_str = record_json_str.strip()
        record_json_str = record_json_str.replace('\n', '')

        return record_json_str

    def pretty_print_record(self):
        record_json_str = JsonUtils.as_json_string(
            JsonUtils.make_writable(self.record)
        )

        print(record_json_str)

    @staticmethod
    def keys() -> list[str]:
        return [
            ORFaqsProteinDiscoveryRecord.ACCESSION_NUMBER_KEY,
            ORFaqsProteinDiscoveryRecord.READING_FRAME_KEY,
            ORFaqsProteinDiscoveryRecord.RNA_SEQUENCE_POSITION_KEY,
            ORFaqsProteinDiscoveryRecord.PROTEIN_KEY,
            ORFaqsProteinDiscoveryRecord.PROTEIN_LENGTH_KEY,
        ]


class ORFaqsProteinDiscoveryUtils:
    """ORFaqsProteinDiscoveryUtils"""

    _MULTITHREADING_SEQUENCE_LENGTH_THRESHOLD = 1000
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
            ORFaqsProteinDiscoveryUtils.DATAFRAME_INDEX_KEY
        ] + ORFaqsProteinDiscoveryRecord.keys()

    @staticmethod
    def unreferenced_accession_number() -> str:
        "unknown_reference"

    @staticmethod
    def _result_file_name(
        accession_number: str = None,
        include_date_time_stamp: bool = False,
        reading_frame: RNAReadingFrame = None,
        custom_tag: str = None,
    ) -> str:
        timestamp_str = datetime.now().strftime('%Y%m%d-%H%M%S')
        file_name = 'discovered-proteins'
        if include_date_time_stamp:
            file_name = f'{timestamp_str}-{file_name}'

        if accession_number is not None:
            file_name = f'{accession_number}-{file_name}'

        if reading_frame is not None:
            file_name = f'{file_name}-reading-frame-{reading_frame.value}'

        if custom_tag is not None:
            file_name = f'{file_name}-{custom_tag}'

        return file_name

    @staticmethod
    def _intermediate_result_file_name(
        reading_frame: RNAReadingFrame,
        custom_tag: str,
    ) -> str:
        file_name = ORFaqsProteinDiscoveryUtils._result_file_name(
            include_date_time_stamp=True,
            reading_frame=reading_frame,
            custom_tag=custom_tag,
        )
        return f'{file_name}.txt'

    @staticmethod
    def _exported_reading_frame_result_file_name(
        access_number: str,
        reading_frame: RNAReadingFrame,
    ) -> str:
        return ORFaqsProteinDiscoveryUtils._result_file_name(
            accession_number=access_number, reading_frame=reading_frame
        )

    @staticmethod
    def exported_file_name(access_number: str = None) -> str:
        return ORFaqsProteinDiscoveryUtils._result_file_name(access_number)

    @staticmethod
    def validate_discovered_proteins_file(file_path: (str | os.PathLike)):
        results_dataframe = PandasUtils.read_file_as_dataframe(file_path)
        expected_columns = ORFaqsProteinDiscoveryRecord.keys()
        for expected_column in expected_columns:
            if expected_column not in results_dataframe.columns:
                message = (
                    f'[ERROR] Expected data column {expected_column} NOT '
                    f'found. Expected columns are: {expected_columns}.'
                )
                _logger.error(message)
                raise ValueError(message)

    @staticmethod
    def is_valid_discovered_proteins_file(
        file_path: (str | os.PathLike),
    ) -> bool:
        try:
            ORFaqsProteinDiscoveryUtils.validate_discovered_proteins_file(
                file_path
            )
        except ValueError:
            return False

        return True

    @staticmethod
    def _read_intermediate_records_file(
        file_path: str | pathlib.Path,
    ) -> list[dict[str, any]]:
        records_list: list[dict[str, any]] = []
        try:
            with open(file_path, 'r', encoding='utf-8') as i_file:
                for line in i_file:
                    records_list.append(JsonUtils.read_json(line.strip()))
        except FileNotFoundError as e:
            message = (
                '[Warning] The file provided could not be found while '
                'running _read_intermediate_records_file. '
                'Continuing...\n'
                '(debug) -> \n'
                '\tThis may be the result of an ERROR, or that an '
                'empty protein list was returned and the file path was '
                'not properly cleared.\n'
                f'{e}'
            )
            _logger.warning(message)
            print(message)

        return records_list

    @staticmethod
    def _export_reading_frame_results(
        accession_number: str,
        output_directory: (str | pathlib.Path),
        reading_frame: RNAReadingFrame,
        file_paths: list[str | pathlib.Path],
        export_format: _ExportFormatOptions,
    ) -> pathlib.Path:
        #######################################################################
        # Read all stored records.
        records_list: list[dict[str, any]] = []
        try:
            for file_path in file_paths:
                records_list += ORFaqsProteinDiscoveryUtils._read_intermediate_records_file(
                    file_path
                )
        except FileNotFoundError as e:
            message = (
                '[ERROR] An input file could not be found while '
                'consolidating results for export.\n'
                f'(debug) ->\n'
                f'\tfile_path: {file_path}'
            )
            _logger.error(message)
            print(message)
            raise FileExistsError from e

        #######################################################################
        # Export the results.

        # Organize the results into a DataFrame object.
        records_dataframe = pd.DataFrame(records_list)

        # Consolidate results for the given reading frame.
        export_file_path = output_directory.joinpath(
            ORFaqsProteinDiscoveryUtils._exported_reading_frame_result_file_name(
                access_number=accession_number,
                reading_frame=reading_frame,
            )
        )
        export_format = str(export_format).lower()
        export_file_path = PandasUtils.export_dataframe(
            file_path=export_file_path,
            dataframe=records_dataframe,
            export_format=export_format,
            index_label=ORFaqsProteinDiscoveryUtils.DATAFRAME_INDEX_KEY,
        )

        #######################################################################
        # Remove intermediate files after all other processing has completed.
        for file_path in file_paths:
            DirectoryUtils.remove_file_path(file_path)

        return export_file_path

    @staticmethod
    def _group_exported_reading_frame_results(
        access_number: str,
        output_directory: (str | pathlib.Path),
        file_paths: list[str | pathlib.Path],
        export_format: _ExportFormatOptions,
    ):
        # Create the dataframe object. Fill it with the
        # expected data column order.
        grouped_results_dataframe: list[pd.DataFrame] = []
        record_keys = ORFaqsProteinDiscoveryRecord.keys()
        for file_path in file_paths:
            if not DirectoryUtils.path_exists(file_path):
                message = (
                    '[WARNING] An expected results file '
                    'could not be found.\n'
                    '(debug) ->\n'
                    f'\tfile_path: {file_path}'
                )
                _logger.warning(message)
                print(message)

            # Read the data from file.
            reading_frame_results = PandasUtils.read_file_as_dataframe(
                file_path
            )
            grouped_results_dataframe.append(
                reading_frame_results[record_keys]
            )
        #######################################################################
        # Export the results.
        final_results_dataframe = pd.concat(
            grouped_results_dataframe,
            axis='index',
            ignore_index=True,
        )

        output_directory = DirectoryUtils.make_path_object(output_directory)
        export_file_path = output_directory.joinpath(
            ORFaqsProteinDiscoveryUtils.exported_file_name(access_number)
        )
        PandasUtils.export_dataframe(
            file_path=export_file_path,
            dataframe=final_results_dataframe,
            export_format=export_format,
            index_label=ORFaqsProteinDiscoveryUtils.DATAFRAME_INDEX_KEY,
        )

    @staticmethod
    def _translate_rna_group(
        accession_number: str,
        reading_frame: RNAReadingFrame,
        rna_sequence: RNASequence,
        start_codon_indices: list[int],
        start_codons: list[Codon],
        stop_codons: list[Codon],
        output_directory: str | pathlib.Path,
        thread_id: int,
        result_queue: multiprocessing.Queue,
        display_progress: bool = True,
    ):
        if len(start_codon_indices) == 0:
            message = '[INFO] EMPTY LIST'
            _logger.info(message)
            print(message)
        protein_list: list[Protein] = []
        # For every start codon found, translate the RNA region
        # beginning at that codon.
        process_list = start_codon_indices
        if display_progress:
            tqdm_description = (
                'Translating RNA from start codons in '
                f'reading frame {reading_frame.value}...'
            )
            if thread_id is not None:
                tqdm_description = f'[Thread {thread_id}] {tqdm_description}'
            process_list = tqdm(
                process_list, desc=tqdm_description, total=len(process_list)
            )
        # Create the intermediate results file.
        results_file_name = (
            ORFaqsProteinDiscoveryUtils._intermediate_result_file_name(
                reading_frame=reading_frame, custom_tag=f'tid{thread_id}'
            )
        )
        results_file_path = DirectoryUtils.make_path_object(
            output_directory
        ).joinpath(results_file_name)
        for start_codon_index in process_list:
            rna_coding_region = rna_sequence[start_codon_index:]
            translated_protein = Ribosome.translate_rna(
                rna_coding_region,
                start_codons=start_codons,
                stop_codons=stop_codons,
            )
            if isinstance(translated_protein, Protein):
                protein_list.append(translated_protein)

            # Record the "1's" based index for the start codon position
            ones_based_index_start_codon_position = start_codon_index + 1
            with open(results_file_path, 'a', encoding='utf-8') as o_file:
                protein_record_json_str = ORFaqsProteinDiscoveryRecord(
                    accession_number=accession_number,
                    reading_frame=reading_frame,
                    rna_sequence_position=ones_based_index_start_codon_position,
                    protein=translated_protein,
                ).condensed_record_json_str
                o_file.write(f'{protein_record_json_str}\n')

        results = (protein_list, results_file_path)
        result_queue.put(results)

    @staticmethod
    def discover_proteins(
        rna_sequence: str | RNASequence,
        accession_number: str = None,
        frame: RNAReadingFrame = None,
        start_codons: list[Codon] = None,
        stop_codons: list[Codon] = None,
        output_directory: str | pathlib.Path = None,
        export_format: _ExportFormatOptions = None,
    ) -> tuple[dict, int]:
        _perf_profiler.start_perf_timer(
            _ProfilingFunctionName.PROTEIN_DISCOVERY_GENOMIC_SEQUENCE
        )
        output_directory = DirectoryUtils.make_path_object(output_directory)
        DirectoryUtils.mkdir_path(output_directory)

        if accession_number is None:
            accession_number = (
                ORFaqsProteinDiscoveryUtils.unreferenced_accession_number()
            )

        if not isinstance(rna_sequence, RNASequence):
            rna_sequence = RNASequence(rna_sequence)

        reading_frames: list[RNAReadingFrame] = []
        if frame is None:
            reading_frames = RibosomeUtils.available_reading_frames()
        else:
            reading_frames = [frame]
        #######################################################################
        # Initialize the protein map and related variables which store the
        # the results of the discovered proteins.
        protein_map = {reading_frame: [] for reading_frame in reading_frames}
        protein_count = 0
        thread_count = 1
        if (
            rna_sequence.sequence_length
            > ORFaqsProteinDiscoveryUtils._MULTITHREADING_SEQUENCE_LENGTH_THRESHOLD
        ):
            thread_count = os.cpu_count()

        exported_reading_frame_file_paths: list[str | pathlib.Path] = []
        for reading_frame in reading_frames:
            (start_index, stop_index) = (
                RibosomeUtils.sequence_start_stop_indices(
                    rna_sequence, reading_frame
                )
            )
            ###################################################################
            # Get all start codon positions for the current reading frame.
            rna_sequence_frame = rna_sequence[start_index:stop_index]
            start_codon_indices = RibosomeUtils.find_start_codons(
                rna_sequence_frame, start_codons
            )
            number_start_codons = len(start_codon_indices)

            if number_start_codons == 0:
                # Skip processing of this reading frame.
                message = (
                    '[INFO] No start codons were found for '
                    f'reading frame {reading_frame.value}.'
                )
                _logger.info(message)
                print(message)
                continue

            ###################################################################
            # Process start codons
            message = (
                f'[INFO] {number_start_codons} start codons found for '
                f'reading frame {reading_frame.value}.'
            )
            _logger.info(message)
            print(message)

            thread_pool: list[multiprocessing.Process] = [None] * thread_count
            return_queue_list = [multiprocessing.Queue()] * thread_count
            group_size = math.ceil(number_start_codons / thread_count)

            for thread_id in range(thread_count):
                start_index = thread_id * group_size
                end_index = start_index + group_size
                start_codon_indices_subset = start_codon_indices[
                    start_index:end_index
                ]

                thread_pool[thread_id] = multiprocessing.Process(
                    target=ORFaqsProteinDiscoveryUtils._translate_rna_group,
                    args=(
                        accession_number,
                        reading_frame,
                        str(rna_sequence_frame),
                        start_codon_indices_subset,
                        start_codons,
                        stop_codons,
                        output_directory,
                        thread_id,
                        return_queue_list[thread_id],
                    ),
                )

                thread_pool[thread_id].start()

            ###################################################################
            # Update the protein map. Because this implementation uses queues,
            # join() is not called. Calling join() with queue objects will
            # cause the application to deadlock.
            result_files_list: list[pathlib.Path] = []
            for return_queue in return_queue_list:
                (protein_list, results_file_path) = return_queue.get()
                # Update the protein map and results file list iff.
                # protein_list is not empty.
                if len(protein_list) > 0:
                    protein_map[reading_frame] += protein_list
                    result_files_list.append(results_file_path)

            protein_count += len(protein_map[reading_frame])
            message = f'[INFO] Protein Count: {protein_count}'
            _logger.info(message)
            print(message)

            ###################################################################
            # Final discovered proteins results export.
            # 1. Create results for the individual reading frames.
            exported_reading_frame_file_path = (
                ORFaqsProteinDiscoveryUtils._export_reading_frame_results(
                    accession_number,
                    output_directory,
                    reading_frame,
                    result_files_list,
                    export_format=export_format,
                )
            )
            exported_reading_frame_file_paths.append(
                exported_reading_frame_file_path
            )

        # 2. Group the results of each reading frame into a single discovered
        # proteins results file.
        if len(exported_reading_frame_file_paths) > 0:
            ORFaqsProteinDiscoveryUtils._group_exported_reading_frame_results(
                accession_number,
                output_directory,
                exported_reading_frame_file_paths,
                export_format=export_format,
            )
        else:
            message = (
                '[INFO] No proteins were found for the given sequence. '
                'No outputs will be generated.'
            )
            _logger.info(message)
            print(message)

        _perf_profiler.stop_perf_timer(
            _ProfilingFunctionName.PROTEIN_DISCOVERY_GENOMIC_SEQUENCE
        )
        return (protein_map, protein_count)

    @staticmethod
    def process_genomic_sequence(
        genomic_sequence: str | GenomicSequence,
        accession_number: str = None,
        export_format: _ExportFormatOptions = None,
        output_directory: str | pathlib.Path = None,
    ):
        if export_format is None:
            export_format = ORFaqsProteinDiscoveryUtils.default_export_format()

        if output_directory is None:
            output_directory = (
                ORFaqsProteinDiscoveryUtils.default_output_directory()
            )

        if isinstance(genomic_sequence, str):
            genomic_sequence = NucleotideUtils.create_sequence(
                genomic_sequence
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

        (protein_map, protein_count) = (
            ORFaqsProteinDiscoveryUtils.discover_proteins(
                rna_sequence,
                accession_number=accession_number,
                output_directory=output_directory,
                export_format=export_format,
            )
        )
        # Write the info file to the output directory

        return (protein_map, protein_count)

    @staticmethod
    def process_fasta_file(
        fasta_file_path: str | pathlib.Path,
        export_format: _ExportFormatOptions = None,
        output_directory: str | pathlib.Path = None,
        display_progress: bool = True,
    ):
        if export_format is None:
            export_format = ORFaqsProteinDiscoveryUtils.default_export_format()

        if output_directory is None:
            output_directory = (
                ORFaqsProteinDiscoveryUtils.default_output_directory()
            )

        fasta_sequences = FASTAUtils.parse_file(fasta_file_path)
        if output_directory is not None:
            output_directory = DirectoryUtils.make_path_object(
                output_directory
            )
        all_protein_maps = []
        protein_counts = {}
        total_protein_count = 0
        process_list = fasta_sequences
        if display_progress:
            tqdm_description = 'Processing genomic sequences...'
            process_list = tqdm(
                process_list, desc=tqdm_description, total=len(process_list)
            )

        for fasta_sequence in process_list:
            print('-------------------------------------------------------')
            print(f'Processing Sequence: {fasta_sequence.sequence_name}')
            current_sequence_output_directory = output_directory.joinpath(
                f'{fasta_sequence.accession_number}'
            )
            (protein_map, protein_count) = (
                ORFaqsProteinDiscoveryUtils.process_genomic_sequence(
                    fasta_sequence.sequence,
                    accession_number=fasta_sequence.accession_number,
                    output_directory=current_sequence_output_directory,
                    export_format=export_format,
                )
            )
            all_protein_maps.append(protein_map)
            protein_counts[fasta_sequence.accession_number] = protein_count
            total_protein_count += protein_count
            print('\n')

        print(f'Total number of proteins discovered: {total_protein_count}')
