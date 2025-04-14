import enum
import logging
import multiprocessing
import os

from datetime import datetime
from enum import Enum
from tqdm import tqdm

from orfaqs.apps.common.cliutils import CliUtil
from orfaqs.apps.common.orfaqsapp import OrfaqsApp

from orfaqs.lib.core.codons import Codon
from orfaqs.lib.core.enzymes import RNAPolymerase
from orfaqs.lib.utils.fastautils import FASTAUtils
from orfaqs.lib.core.nucleotides import (
    DNASequence,
    GenomicSequence,
    NucleotideUtils,
    RNASequence
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


class _ProfilingFunctionName(Enum):
    PROTEIN_DISCOVERY_GENOMIC_SEQUENCE = enum.auto()
    TRANSLATE_RNA_GROUP = enum.auto()


class _ProteinDiscoveryRecord:
    '''_ProteinDiscoveryRecord'''
    _READING_FRAME_KEY = 'reading_frame'
    _RNA_SEQUENCE_POSITION = 'rna_sequence_position'
    _PROTEIN = 'protein'
    _PROTEIN_LENGTH = 'protein_length'

    def __init__(self,
                 reading_frame: RNAReadingFrame,
                 rna_sequence_position: int,
                 protein: Protein):
        self._reading_frame = reading_frame
        self._rna_sequence_position = rna_sequence_position
        self._protein = protein

    @property
    def record(self) -> dict:
        return {
            _ProteinDiscoveryRecord._READING_FRAME_KEY: self._reading_frame.value,
            _ProteinDiscoveryRecord._RNA_SEQUENCE_POSITION: (
                self._rna_sequence_position),
            _ProteinDiscoveryRecord._PROTEIN: self._protein,
            _ProteinDiscoveryRecord._PROTEIN_LENGTH: (
                self._protein.number_amino_acids)
        }

    @property
    def condensed_record_json_str(self) -> str:
        record_json_str = JsonUtils.as_json_string(
            JsonUtils.make_writable(self.record),
            indent=None
        )
        record_json_str = record_json_str.strip()
        record_json_str = record_json_str.replace('\n', '')

        return record_json_str

    def pretty_print_record(self):
        record_json_str = JsonUtils.as_json_string(
            JsonUtils.make_writable(self.record)
        )

        print(record_json_str)


class OrfaqsProteinDiscovery(OrfaqsApp):
    '''OrfaqsProteinDiscovery'''
    _MULTITHREADING_SEQUENCE_LENGTH_THRESHOLD = 1000

    @staticmethod
    def _translate_rna_group(reading_frame: RNAReadingFrame,
                             rna_sequence: RNASequence,
                             start_codon_indices: list[int],
                             start_codons: list[Codon],
                             stop_codons: list[Codon],
                             log_file_path: str | os.PathLike,
                             thread_index: int,
                             result_queue: multiprocessing.Queue):

        protein_list: list[Protein] = []
        log_file_path = DirectoryUtils.make_path_object(log_file_path)
        # For every start codon found, translate the RNA region
        # beginning at that codon.
        process_description = ('Translating RNA from start codons in '
                               f'reading frame {reading_frame.value}...')
        if thread_index is not None:
            process_description = (f'[Thread {thread_index}] '
                                   f'{process_description}')
            parent_directory = log_file_path.parent
            log_file_name = DirectoryUtils.remove_extension(log_file_path.name)
            log_file_extension = log_file_path.suffix
            log_file_path = parent_directory.joinpath(
                (f'{log_file_name}-tid{thread_index}')
            ).with_suffix(log_file_extension)

        for start_codon_index in tqdm(start_codon_indices,
                                      desc=process_description):
            rna_coding_region = rna_sequence[start_codon_index:]
            translated_protein = Ribosome.translate_rna(
                rna_coding_region,
                start_codons=start_codons,
                stop_codons=stop_codons
            )
            if isinstance(translated_protein, Protein):
                protein_list.append(translated_protein)

            # Record the "1's" based index for the start codon position
            ones_based_index_start_codon_position = start_codon_index + 1
            with open(log_file_path, 'a', encoding='utf-8') as i_file:
                protein_record_json_str = _ProteinDiscoveryRecord(
                    reading_frame,
                    ones_based_index_start_codon_position,
                    translated_protein
                ).condensed_record_json_str
                i_file.write(f'{protein_record_json_str}\n')

        result_queue.put(protein_list)

    @staticmethod
    def discover_proteins(
            rna_sequence: str | RNASequence,
            frame: RNAReadingFrame = None,
            start_codons: list[Codon] = None,
            stop_codons: list[Codon] = None,
            log_directory_path: str | os.PathLike = None) -> tuple[dict, int]:

        _perf_profiler.start_perf_timer(
            _ProfilingFunctionName.PROTEIN_DISCOVERY_GENOMIC_SEQUENCE
        )
        if log_directory_path is None:
            log_directory_path = '.log-path'

        log_directory_path = DirectoryUtils.make_path_object(
            log_directory_path
        )
        DirectoryUtils.mkdir_path(log_directory_path)
        timestamp_str = datetime.now().strftime('%Y%m%d-%H%M%S')
        file_name = f'{timestamp_str}-discovered-proteins.txt'
        log_file_path = log_directory_path.joinpath(file_name)

        if not isinstance(rna_sequence, RNASequence):
            rna_sequence = RNASequence(rna_sequence)

        reading_frames: list[RNAReadingFrame] = []
        if frame is None:
            reading_frames = RibosomeUtils.available_reading_frames()
        else:
            reading_frames = [frame]

        protein_map = {reading_frame: [] for reading_frame in reading_frames}
        protein_count = 0
        thread_count = 1
        if (rna_sequence.sequence_length >
                OrfaqsProteinDiscovery._MULTITHREADING_SEQUENCE_LENGTH_THRESHOLD):
            thread_count = os.cpu_count()

        for reading_frame in reading_frames:
            (start_index,
             stop_index) = RibosomeUtils.sequence_start_stop_indices(
                rna_sequence,
                reading_frame
            )
            rna_sequence_frame = rna_sequence[start_index:stop_index]
            # Get all start codon positions for the given reading frame.
            start_codon_indices = RibosomeUtils.find_start_codons(
                rna_sequence_frame,
                start_codons
            )
            number_start_codons = len(start_codon_indices)
            thread_pool: list[multiprocessing.Process] = [None] * thread_count
            return_queue_list = [multiprocessing.Queue()] * thread_count
            group_size = int(number_start_codons / thread_count)

            # Process start codons
            message = (f'{number_start_codons} start codons found for '
                       f'reading frame {reading_frame.value}.')
            _logger.info(message)
            print(message)
            for thread_index in range(thread_count):
                start_index = thread_index*group_size
                end_index = start_index + group_size
                start_codon_indices_subset = (
                    start_codon_indices[start_index:end_index]
                )
                thread_pool[thread_index] = multiprocessing.Process(
                    target=OrfaqsProteinDiscovery._translate_rna_group,
                    args=(
                        reading_frame,
                        str(rna_sequence_frame),
                        start_codon_indices_subset,
                        start_codons,
                        stop_codons,
                        log_file_path,
                        thread_index,
                        return_queue_list[thread_index]
                    )
                )

                thread_pool[thread_index].start()

            # Update the protein map. Because this implementation uses queues,
            # join() is not called. Calling join() with queue objects will
            # cause the application to deadlock.
            for return_queue in return_queue_list:
                protein_map[reading_frame] += return_queue.get()

            protein_count += len(protein_map[reading_frame])
            print(f'Protein Count: {protein_count}')

        _perf_profiler.stop_perf_timer(
            _ProfilingFunctionName.PROTEIN_DISCOVERY_GENOMIC_SEQUENCE
        )
        return (protein_map, protein_count)

    @staticmethod
    def _process_genomic_sequence(
            genomic_sequence: str | GenomicSequence,
            log_directory_path: str | os.PathLike = None):
        if isinstance(genomic_sequence, str):
            genomic_sequence = NucleotideUtils.create_sequence(
                genomic_sequence
            )

        rna_sequence: RNASequence = None
        if isinstance(genomic_sequence, DNASequence):
            rna_sequence = RNAPolymerase.transcribe(genomic_sequence)
        elif isinstance(genomic_sequence, RNASequence):
            rna_sequence = genomic_sequence

        (protein_map,
         protein_count) = OrfaqsProteinDiscovery.discover_proteins(
             rna_sequence,
             log_directory_path=log_directory_path
        )
        return (protein_map, protein_count)

    @staticmethod
    def _process_fasta_file(fasta_file_path: str | os.PathLike,
                            output_directory: str | os.PathLike = None):
        fasta_sequences = FASTAUtils.parse_file(fasta_file_path)
        if output_directory is not None:
            output_directory = DirectoryUtils.make_path_object(
                output_directory
            )
        all_protein_maps = []
        protein_counts = {}
        total_protein_count = 0
        for (sequence_id, fasta_sequence) in enumerate(fasta_sequences):
            print('-------------------------------------------------------')
            print(f'Processing Sequence: {fasta_sequence.sequence_name}')
            log_directory = output_directory.joinpath(
                f'sequence-{sequence_id}'
            )
            (protein_map,
             protein_count) = OrfaqsProteinDiscovery.process_genomic_sequence(
                 fasta_sequence.sequence,
                 log_directory_path=log_directory
            )
            all_protein_maps.append(protein_map)
            protein_counts[fasta_sequence.sequence_name] = protein_count
            total_protein_count += protein_count
            print('\n')

        print(f'Total number of proteins discovered: {total_protein_count}')

    @staticmethod
    def program_name():
        return 'ORFAQs Protein Discovery'

    @staticmethod
    def cli():
        try:
            arg_descriptor_list = [
                CliUtil.create_new_arg_descriptor(
                    ('-j', '--job_id'),
                    arg_help=('A string used to uniquely identify a set of '
                              'results generated by the tool. All results '
                              'with the same job_id share the same directory '
                              'sub-folder.'),
                    arg_type=str,
                    default=None),
                CliUtil.create_new_arg_descriptor(
                    ('-i', '--input_sequence'),
                    arg_help=('A sequence file or genomic sequence string. '
                              'Currently, only FASTA files are supported.'),
                    arg_type=str),
                CliUtil.create_new_arg_descriptor(
                    ('-o', '--output_directory'),
                    arg_help=('The output directory for exported files '
                              'and results. If the path does not exist, one '
                              'will be created. The default directory is: '
                              f'{OrfaqsProteinDiscovery.default_output_directory()}'),
                    arg_type=str,
                    default=OrfaqsProteinDiscovery.default_output_directory())
            ]

            cli_arg_parser = CliUtil.create_arg_parser(
                arg_descriptor_list,
                program_name=OrfaqsProteinDiscovery.program_name(),
                description=('Provide block level collections of images.'),
                epilog='')

            ui_args = CliUtil.parse_args(cli_arg_parser)
            if not ui_args['input_sequence']:
                info_message = ('No input directory or sequence was '
                                'specified. An input must be given for the '
                                'program to generate its outputs.')
                print(info_message)
                cli_arg_parser.print_help()
            else:
                OrfaqsProteinDiscovery.run(**ui_args)

        except SystemExit as error:
            if error.args[0] is True:
                raise

    @staticmethod
    def _run(input_sequence: str | os.PathLike,
             output_directory: str | os.PathLike,
             job_id: str = None,
             **kwargs):

        if ((input_sequence is None) and
            (not DirectoryUtils.path_exists(input_sequence)) and
            (not NucleotideUtils.is_dna_sequence(input_sequence)) and
                (not NucleotideUtils.is_rna_sequence(input_sequence))):
            message = ('[Invalid Argument]: A valid input object must be '
                       f'specified when calling {__class__.__qualname__}')
            _logger.error(message)
            raise ValueError(message)

        if output_directory is None:
            output_directory = OrfaqsProteinDiscovery.default_output_directory()

        # Create the local output directory and output file path
        output_directory = DirectoryUtils.make_path_object(output_directory)
        if isinstance(job_id, str):
            output_directory = output_directory.joinpath(job_id)

        if DirectoryUtils.is_file(input_sequence):
            # Process as a FASTA file
            OrfaqsProteinDiscovery._process_fasta_file(
                input_sequence,
                output_directory
            )
        else:
            OrfaqsProteinDiscovery._process_genomic_sequence(
                input_sequence,
                output_directory
            )


if __name__ == '__main__':
    OrfaqsProteinDiscovery.cli()
