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
import typing

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
from orfaqs.lib.utils.perfutils import PerfProfiler

_logger = logging.getLogger(__name__)
_perf_profiler = PerfProfiler(__name__)


class _ProfilingFunctionName(enum.Enum):
    PROTEIN_DISCOVERY_GENOMIC_SEQUENCE = enum.auto()
    TRANSLATE_RNA_GROUP = enum.auto()


class _ExportFormat(enum.Enum):
    CSV = 'csv'
    JSON = 'json'
    XLSX = 'xlsx'


_AVAILABLE_EXPORT_FORMATS = [format.value for format in _ExportFormat]
_EXPORT_FORMAT_OPTIONS = typing.Literal[*tuple(_AVAILABLE_EXPORT_FORMATS)]


class ORFaqsProteinDiscoveryRecord:
    """ORFaqsProteinDiscoveryRecord"""

    PROTEIN_KEY = 'protein'
    PROTEIN_LENGTH_KEY = 'protein_length'
    READING_FRAME_KEY = 'reading_frame'
    RNA_SEQUENCE_POSITION_KEY = 'rna_sequence_position'

    def __init__(
        self,
        reading_frame: RNAReadingFrame,
        rna_sequence_position: int,
        protein: Protein,
    ):
        self._reading_frame = reading_frame
        self._rna_sequence_position = rna_sequence_position
        self._protein = protein

    @property
    def record(self) -> dict[str, any]:
        record_map: dict[str, any] = {}
        for record_key in ORFaqsProteinDiscoveryRecord.keys():
            if ORFaqsProteinDiscoveryRecord.READING_FRAME_KEY == record_key:
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
        return _ExportFormat.CSV.value

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
    def read_results_file_as_dataframe(
        results_file: (str | os.PathLike),
    ) -> pd.DataFrame:
        results_dataframe: pd.DataFrame = None
        exported_file_type = DirectoryUtils.make_path_object(
            results_file
        ).suffix.lower()
        if _ExportFormat.CSV.value in exported_file_type:
            results_dataframe = pd.read_csv(results_file, index_col=0)
        elif _ExportFormat.JSON.value in exported_file_type:
            results_dataframe = pd.read_json(results_file)
        elif _ExportFormat.XLSX.value in exported_file_type:
            results_dataframe = pd.read_excel(results_file, index_col=0)

        return results_dataframe

    @staticmethod
    def _result_file_name(
        include_date_time_stamp: bool = False,
        reading_frame: RNAReadingFrame = None,
        custom_tag: str = None,
    ) -> pathlib.Path:
        timestamp_str = datetime.now().strftime('%Y%m%d-%H%M%S')
        file_name = 'discovered-proteins'
        if include_date_time_stamp:
            file_name = f'{timestamp_str}-{file_name}'

        if reading_frame is not None:
            file_name = f'{file_name}-reading-frame-{reading_frame.value}'

        if custom_tag is not None:
            file_name = f'{file_name}-{custom_tag}'

        return DirectoryUtils.make_path_object(file_name).with_suffix('.txt')

    @staticmethod
    def _intermediate_result_file_name(
        reading_frame: RNAReadingFrame, custom_tag: str
    ) -> pathlib.Path:
        return ORFaqsProteinDiscoveryUtils._result_file_name(
            include_date_time_stamp=True,
            reading_frame=reading_frame,
            custom_tag=custom_tag,
        )

    @staticmethod
    def _exported_result_file_name(
        reading_frame: RNAReadingFrame,
    ) -> pathlib.Path:
        return ORFaqsProteinDiscoveryUtils._result_file_name(
            include_date_time_stamp=False,
            reading_frame=reading_frame,
            custom_tag=None,
        )

    @staticmethod
    def _final_result_file_name(custom_tag: str = None) -> pathlib.Path:
        return ORFaqsProteinDiscoveryUtils._result_file_name(
            custom_tag=custom_tag
        )

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
    def _export_intermediate_results(
        output_file_path: (str | pathlib.Path),
        file_paths: list[str | pathlib.Path],
        export_format: _EXPORT_FORMAT_OPTIONS,
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
        # Organize the results into a DataFrame object.
        records_dataframe = pd.DataFrame(records_list)
        #######################################################################
        # Export the results.
        export_format = str(export_format).lower()
        export_file_path: pathlib.Path = None
        if _ExportFormat.CSV.value in export_format:
            export_file_path = output_file_path.with_suffix(
                f'.{_ExportFormat.CSV.value}'
            )
            records_dataframe.to_csv(
                export_file_path,
                index_label=ORFaqsProteinDiscoveryUtils.DATAFRAME_INDEX_KEY,
            )
        elif _ExportFormat.JSON.value in export_format:
            export_file_path = output_file_path.with_suffix(
                f'.{_ExportFormat.JSON.value}'
            )
            records_dataframe.to_json(export_file_path)
        elif _ExportFormat.XLSX.value in export_format:
            export_file_path = output_file_path.with_suffix(
                f'.{_ExportFormat.XLSX.value}'
            )
            records_dataframe.to_excel(
                export_file_path,
                index_label=ORFaqsProteinDiscoveryUtils.DATAFRAME_INDEX_KEY,
            )

        #######################################################################
        # Remove intermediate files after all other
        # processing has completed successfully.
        for file_path in file_paths:
            DirectoryUtils.remove_file_path(file_path)

        return export_file_path

    @staticmethod
    def _group_exported_reading_frame_results(
        output_file_path: (str | pathlib.Path),
        file_paths: list[str | pathlib.Path],
        export_file_type: _ExportFormat,
    ):
        # Create the dataframe object. Fill it with the
        # expected data column order.
        grouped_results_dataframe: pd.DataFrame = pd.DataFrame()
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

            reading_frame_results: pd.DataFrame = None
            if _ExportFormat.CSV.value in export_file_type:
                reading_frame_results = pd.read_csv(file_path, index_col=0)
            elif _ExportFormat.JSON.value in export_file_type:
                reading_frame_results = pd.read_json(file_path)
            elif _ExportFormat.XLSX.value in export_file_type:
                reading_frame_results = pd.read_excel(file_path, index_col=0)
            reading_frame_results = reading_frame_results.reset_index()
            reading_frame_results = reading_frame_results[record_keys]
            # reading_frame_results = reading_frame_results[expected_keys]
            grouped_results_dataframe = pd.concat(
                [grouped_results_dataframe, reading_frame_results],
                axis=0,
                ignore_index=True,
            )
        #######################################################################
        # Export the results.
        if _ExportFormat.CSV.value in export_file_type:
            grouped_results_dataframe.to_csv(
                output_file_path.with_suffix(f'.{_ExportFormat.CSV.value}'),
                index_label=ORFaqsProteinDiscoveryUtils.DATAFRAME_INDEX_KEY,
            )
        elif _ExportFormat.JSON.value in export_file_type:
            grouped_results_dataframe.to_json(
                output_file_path.with_suffix(f'.{_ExportFormat.JSON.value}')
            )
        elif _ExportFormat.XLSX.value in export_file_type:
            grouped_results_dataframe.to_excel(
                output_file_path.with_suffix(f'.{_ExportFormat.XLSX.value}'),
                index_label=ORFaqsProteinDiscoveryUtils.DATAFRAME_INDEX_KEY,
            )

    @staticmethod
    def translate_rna_group(
        reading_frame: RNAReadingFrame,
        rna_sequence: RNASequence,
        start_codon_indices: list[int],
        start_codons: list[Codon],
        stop_codons: list[Codon],
        output_directory: str | pathlib.Path,
        thread_index: int,
        result_queue: multiprocessing.Queue,
    ):
        if len(start_codon_indices) == 0:
            message = '[INFO] EMPTY LIST'
            _logger.info(message)
            print(message)
        protein_list: list[Protein] = []
        # For every start codon found, translate the RNA region
        # beginning at that codon.
        process_list = tqdm(
            start_codon_indices,
            desc=(
                f'[Thread {thread_index}] '
                'Translating RNA from start codons in '
                f'reading frame {reading_frame.value}...'
            ),
        )
        # Create the intermediate results file.
        results_file_name = (
            ORFaqsProteinDiscoveryUtils._intermediate_result_file_name(
                reading_frame=reading_frame, custom_tag=f'tid{thread_index}'
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
                    reading_frame,
                    ones_based_index_start_codon_position,
                    translated_protein,
                ).condensed_record_json_str
                o_file.write(f'{protein_record_json_str}\n')

        results = (protein_list, results_file_path)
        result_queue.put(results)

    @staticmethod
    def discover_proteins(
        rna_sequence: str | RNASequence,
        frame: RNAReadingFrame = None,
        start_codons: list[Codon] = None,
        stop_codons: list[Codon] = None,
        output_directory: str | pathlib.Path = None,
        export_format: _ExportFormat = None,
    ) -> tuple[dict, int]:
        _perf_profiler.start_perf_timer(
            _ProfilingFunctionName.PROTEIN_DISCOVERY_GENOMIC_SEQUENCE
        )
        output_directory = DirectoryUtils.make_path_object(output_directory)
        DirectoryUtils.mkdir_path(output_directory)

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

            for thread_index in range(thread_count):
                start_index = thread_index * group_size
                end_index = start_index + group_size
                start_codon_indices_subset = start_codon_indices[
                    start_index:end_index
                ]

                thread_pool[thread_index] = multiprocessing.Process(
                    target=ORFaqsProteinDiscoveryUtils.translate_rna_group,
                    args=(
                        reading_frame,
                        str(rna_sequence_frame),
                        start_codon_indices_subset,
                        start_codons,
                        stop_codons,
                        output_directory,
                        thread_index,
                        return_queue_list[thread_index],
                    ),
                )

                thread_pool[thread_index].start()

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

            # Consolidate results for the given reading frame.
            reading_frame_results_file_name = (
                ORFaqsProteinDiscoveryUtils._exported_result_file_name(
                    reading_frame=reading_frame
                )
            )
            final_results_file_path = output_directory.joinpath(
                reading_frame_results_file_name
            )
            exported_reading_frame_file_path = (
                ORFaqsProteinDiscoveryUtils._export_intermediate_results(
                    final_results_file_path,
                    result_files_list,
                    export_format=export_format,
                )
            )
            exported_reading_frame_file_paths.append(
                exported_reading_frame_file_path
            )

        final_output_file_name = (
            ORFaqsProteinDiscoveryUtils._final_result_file_name()
        )
        final_output_file_path = output_directory.joinpath(
            final_output_file_name
        )

        if len(exported_reading_frame_file_paths) > 0:
            ORFaqsProteinDiscoveryUtils._group_exported_reading_frame_results(
                final_output_file_path,
                exported_reading_frame_file_paths,
                export_file_type=export_format,
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
        export_format: _ExportFormat = None,
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
                output_directory=output_directory,
                export_format=export_format,
            )
        )
        return (protein_map, protein_count)

    @staticmethod
    def process_fasta_file(
        fasta_file_path: str | pathlib.Path,
        export_format: _EXPORT_FORMAT_OPTIONS = None,
        output_directory: str | pathlib.Path = None,
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
        for fasta_sequence in fasta_sequences:
            print('-------------------------------------------------------')
            print(f'Processing Sequence: {fasta_sequence.sequence_name}')
            current_sequence_output_directory = output_directory.joinpath(
                f'sequence-{fasta_sequence.sequence_id}'
            )
            (protein_map, protein_count) = (
                ORFaqsProteinDiscoveryUtils.process_genomic_sequence(
                    fasta_sequence.sequence,
                    output_directory=current_sequence_output_directory,
                    export_format=export_format,
                )
            )
            all_protein_maps.append(protein_map)
            protein_counts[fasta_sequence.sequence_name] = protein_count
            total_protein_count += protein_count
            print('\n')

        print(f'Total number of proteins discovered: {total_protein_count}')
