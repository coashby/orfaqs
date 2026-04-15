"""
ORFaqs Protein Query common app classes, resources, and utility functions.
"""

import logging
import os
import pandas as pd
import pathlib
import typing

from orfaqs.apps.common.orfaqsproteindiscovery import (
    RNAReadingFrame,
)
from orfaqs.apps.common.orfaqsrecords import (
    ORFaqsDiscoveredProteinRecord,
    ORFaqsDiscoveredProteinRecordKeys,
    ORFaqsProteinRecord,
    ORFaqsRecordUtils,
)

from orfaqs.lib.core.nucleotides import NucleotideUtils
from orfaqs.lib.utils.databaseutils import (
    BaseTable,
    PostgresDatabaseUtils,
    sqlalchemy,
    SqlAlchemyUtils,
)
from orfaqs.lib.utils.directoryutils import DirectoryUtils
from orfaqs.lib.utils.fastautils import (
    FASTASequenceType,
    FASTAUtils,
)
from orfaqs.lib.utils.pandasutils import (
    DataFrameExportFormat,
    PandasUtils,
)

_logger = logging.getLogger(__name__)

_ExportFormatOptions = typing.Literal[
    DataFrameExportFormat.CSV,
    DataFrameExportFormat.XLSX,
]
_AVAILABLE_EXPORT_FORMATS: list[str] = [
    format for format in _ExportFormatOptions.__args__
]

_ORFAQS_PROTEIN_QUERY_DATABASE = 'orfaqs_protein_query'
# Global database connection options enable user defined options to persist
# during the entirety of a session.
_database_connection_options = (
    PostgresDatabaseUtils.default_database_connection_options()
)
_database_connection_options.database = _ORFAQS_PROTEIN_QUERY_DATABASE


class ORFaqsDiscoveredProteinTableSchema(ORFaqsDiscoveredProteinRecordKeys):
    UID_KEY = 'uid'

    @staticmethod
    def columns() -> list[str]:
        return [
            ORFaqsDiscoveredProteinTableSchema.UID_KEY,
            ORFaqsDiscoveredProteinTableSchema.SOURCE_UID_KEY,
            ORFaqsDiscoveredProteinTableSchema.STRAND_TYPE_KEY,
            ORFaqsDiscoveredProteinTableSchema.READING_FRAME_KEY,
            ORFaqsDiscoveredProteinTableSchema.DNA_SEQUENCE_POSITION_KEY,
            ORFaqsDiscoveredProteinTableSchema.PROTEIN_KEY,
            ORFaqsDiscoveredProteinTableSchema.PROTEIN_LENGTH_KEY,
        ]


class _ORFaqsProteinTableUtils:
    @staticmethod
    def uid_comment() -> str:
        return 'The unique id of the entry.'

    @staticmethod
    def source_uid_comment() -> str:
        return (
            'The identifier assigned to the sequence from '
            'which the protein was translated.'
        )

    @staticmethod
    def strand_type_comment() -> str:
        strand_types_str = NucleotideUtils.available_strand_types_str()
        return (
            f'The genomic sequence strand {strand_types_str} used '
            'during translation.'
        )

    @staticmethod
    def reading_frame_comment() -> str:
        return (
            'The reading frame (1, 2, or 3) from which the '
            'protein was transcribed.'
        )

    @staticmethod
    def rna_sequence_position_comment() -> str:
        return (
            'The position of the start codon in the RNA sequence '
            'from which the protein was transcribed.'
        )

    @staticmethod
    def protein_comment() -> str:
        return (
            'The protein represented by its single-letter '
            'amino acid abbreviations.'
        )

    @staticmethod
    def protein_length_comment() -> str:
        return 'The number of amino acids in the protein sequence.'


class _ORFaqsReferenceProteinsTableFactory:
    """_ORFaqsReferenceProteinsTableFactory"""

    TABLE_NAME = 'reference_proteins'

    @staticmethod
    def define_table():
        """
        Defines the discovered_proteins table when called iff the table
        definition does NOT exists in the base class.
        """

        if (
            _ORFaqsReferenceProteinsTableFactory.TABLE_NAME
            in BaseTable.metadata.tables
        ):
            return

        class _ORFaqsReferenceProteinsTable(BaseTable):
            """_ORFaqsReferenceProteinsTable"""

            __tablename__ = _ORFaqsReferenceProteinsTableFactory.TABLE_NAME
            UID_MAX_CHAR_LENGTH = 64

            @staticmethod
            def _reading_frame_check_constraint() -> str:
                reading_frame_values = [
                    reading_frame.value for reading_frame in RNAReadingFrame
                ]
                return f'reading_frame IN {tuple(reading_frame_values)}'

            uid = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.UID_KEY,
                sqlalchemy.String(UID_MAX_CHAR_LENGTH),
                primary_key=True,
                comment=_ORFaqsProteinTableUtils.uid_comment(),
            )

            protein = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.PROTEIN_KEY,
                sqlalchemy.String,
                nullable=False,
                comment=_ORFaqsProteinTableUtils.protein_comment(),
            )

            protein_length = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.PROTEIN_LENGTH_KEY,
                sqlalchemy.Integer,
                nullable=False,
                comment=_ORFaqsProteinTableUtils.protein_length_comment(),
            )


class _ORFaqsDiscoveredProteinsTableFactory:
    """_ORFaqsDiscoveredProteinsTableUtil"""

    TABLE_NAME = 'discovered_proteins'

    @staticmethod
    def define_table():
        """
        Defines the discovered_proteins table when called iff the table
        definition does NOT exists in the base class.
        """
        if (
            _ORFaqsDiscoveredProteinsTableFactory.TABLE_NAME
            in BaseTable.metadata.tables
        ):
            return

        class ORFaqsDiscoveredProteinsTable(BaseTable):
            """ORFaqsDiscoveredProteinTable"""

            __tablename__ = _ORFaqsDiscoveredProteinsTableFactory.TABLE_NAME
            UID_MAX_CHAR_LENGTH = 64
            SOURCE_UID_MAX_CHAR_LENGTH = 64
            STRAND_TYPE_MAX_CHAR_LENGTH = 32

            @staticmethod
            def _reading_frame_check_constraint() -> str:
                reading_frame_values = [
                    reading_frame.value for reading_frame in RNAReadingFrame
                ]
                return f'reading_frame IN {tuple(reading_frame_values)}'

            uid = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.UID_KEY,
                sqlalchemy.String(UID_MAX_CHAR_LENGTH),
                primary_key=True,
                comment=_ORFaqsProteinTableUtils.uid_comment(),
            )

            source_uid = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.SOURCE_UID_KEY,
                sqlalchemy.String(SOURCE_UID_MAX_CHAR_LENGTH),
                comment=_ORFaqsProteinTableUtils.source_uid_comment(),
            )

            strand_type = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.STRAND_TYPE_KEY,
                sqlalchemy.String(STRAND_TYPE_MAX_CHAR_LENGTH),
                nullable=False,
                comment=_ORFaqsProteinTableUtils.strand_type_comment(),
            )

            reading_frame = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.READING_FRAME_KEY,
                sqlalchemy.Integer,
                nullable=False,
                comment=_ORFaqsProteinTableUtils.reading_frame_comment(),
            )

            rna_sequence_position = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.DNA_SEQUENCE_POSITION_KEY,
                sqlalchemy.Integer,
                sqlalchemy.CheckConstraint(_reading_frame_check_constraint()),
                nullable=False,
                comment=_ORFaqsProteinTableUtils.rna_sequence_position_comment(),
            )

            protein = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.PROTEIN_KEY,
                sqlalchemy.String,
                nullable=False,
                comment=_ORFaqsProteinTableUtils.protein_comment(),
            )

            protein_length = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.PROTEIN_LENGTH_KEY,
                sqlalchemy.Integer,
                nullable=False,
                comment=_ORFaqsProteinTableUtils.protein_length_comment(),
            )


class ORFaqsProteinQueryApi:
    """ORFaqsProteinQueryApi"""

    @staticmethod
    def default_export_format() -> str:
        return DataFrameExportFormat.CSV

    @staticmethod
    def default_output_directory() -> pathlib.Path:
        return DirectoryUtils.make_path_object('./')

    @staticmethod
    def default_export_file_name() -> pathlib.Path:
        return 'orfaqs-protein-query-results'

    @staticmethod
    def default_export_file_path() -> pathlib.Path:
        default_format = ORFaqsProteinQueryApi.default_export_format()
        default_file_name = ORFaqsProteinQueryApi.default_export_file_name()
        return (
            ORFaqsProteinQueryApi.default_output_directory()
            / f'{default_file_name}.{default_format}'
        )

    @staticmethod
    def available_export_formats() -> list[str]:
        return _AVAILABLE_EXPORT_FORMATS

    @staticmethod
    def default_database():
        return _ORFAQS_PROTEIN_QUERY_DATABASE

    @staticmethod
    def set_database(database: str):
        global _database_connection_options
        _database_connection_options.database = database

    @staticmethod
    def set_workspace(workspace: str):
        ORFaqsProteinQueryApi.set_database(workspace)

    @staticmethod
    def table() -> str:
        return _ORFaqsDiscoveredProteinsTableFactory.TABLE_NAME

    @staticmethod
    def configure_database(
        username: str = None,
        password: str = None,
        host: str = None,
        port: str = None,
    ):
        global _database_connection_options
        _database_connection_options.username = username
        _database_connection_options.password = password
        _database_connection_options.host = host
        _database_connection_options.port = port

    @staticmethod
    def _create_orfaqs_discovered_proteins_table():
        _ORFaqsDiscoveredProteinsTableFactory.define_table()
        PostgresDatabaseUtils.create_table(
            BaseTable,
            _database_connection_options,
        )

    @staticmethod
    def _create_orfaqs_reference_proteins_table():
        _ORFaqsReferenceProteinsTableFactory.define_table()
        PostgresDatabaseUtils.create_table(
            BaseTable,
            _database_connection_options,
        )

    @staticmethod
    def _create_uid(
        uid: str,
        strand_type: str,
        reading_frame: int,
        rna_sequence_position: int,
        protein_length: str,
    ) -> str:
        uid_str = ':'.join(
            [
                uid,
                strand_type,
                str(reading_frame),
                str(rna_sequence_position),
                str(protein_length),
            ]
        )
        return uid_str

    @staticmethod
    def _add_database_primary_key_column(proteins_dataframe: pd.DataFrame):
        def create_uid(row: pd.Series):
            uid = row[ORFaqsDiscoveredProteinTableSchema.SOURCE_UID_KEY]
            strand_type = row[
                ORFaqsDiscoveredProteinTableSchema.STRAND_TYPE_KEY
            ]
            reading_frame = row[
                ORFaqsDiscoveredProteinTableSchema.READING_FRAME_KEY
            ]
            rna_sequence_position = row[
                ORFaqsDiscoveredProteinTableSchema.DNA_SEQUENCE_POSITION_KEY
            ]
            protein_length = row[
                ORFaqsDiscoveredProteinTableSchema.PROTEIN_LENGTH_KEY
            ]
            return ORFaqsProteinQueryApi._create_uid(
                uid=uid,
                strand_type=strand_type,
                reading_frame=reading_frame,
                rna_sequence_position=rna_sequence_position,
                protein_length=protein_length,
            )

        proteins_dataframe[ORFaqsDiscoveredProteinTableSchema.UID_KEY] = (
            proteins_dataframe.apply(create_uid, axis='columns')
        )

    @staticmethod
    def _prepare_dataframe_for_table(proteins_dataframe: pd.DataFrame):
        # Create a primary key
        ORFaqsProteinQueryApi._add_database_primary_key_column(
            proteins_dataframe
        )
        expected_columns = ORFaqsDiscoveredProteinTableSchema.columns()
        drop_columns = proteins_dataframe.columns.difference(expected_columns)
        proteins_dataframe.drop(drop_columns, axis='columns', inplace=True)

    @staticmethod
    def _prepare_dataframe_for_export(proteins_dataframe: pd.DataFrame):
        proteins_dataframe.drop(
            columns=[ORFaqsDiscoveredProteinTableSchema.UID_KEY],
            inplace=True,
            errors='ignore',
        )

    @staticmethod
    def create_workspace(workspace: str):
        ORFaqsProteinQueryApi.set_workspace(workspace)
        PostgresDatabaseUtils.create_database(_database_connection_options)

    @staticmethod
    def remove_workspace(workspace: str):
        PostgresDatabaseUtils.drop_database(workspace)

    @staticmethod
    def load_discovered_proteins(
        workspace: str,
        input_path: (str | os.PathLike),
    ):
        """
        Reads the contents of the file or directory path provided and loads the
        data into the defined database. If the provided path is a directory, it
        must contain the file ORFaqsProteinDiscoveryUtils.exported_file_name()

        Attempts to load data are only made if data checks are passed.
        Otherwise, the process is aborted and no data is available.

        ---------
        Arguments
        ---------
        workspace (str):
            The name of the workspace session to load the discovered proteins.
            Workspaces persist until they are removed by the calling
            application.

        input_path (str | os.PathLike):
            A file or directory path to discovered proteins.
        """
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

        discovered_proteins_files: list[os.PathLike] = []
        for file_path in input_file_paths:
            # Validate the file paths...
            try:
                proteins_dataframe = PandasUtils.read_file_as_dataframe(
                    file_path
                )
                for record_key in ORFaqsDiscoveredProteinRecord.keys():
                    if record_key not in proteins_dataframe.columns:
                        continue
            except ValueError:
                continue
            discovered_proteins_files.append(file_path)

        if len(discovered_proteins_files) == 0:
            message = (
                '[INFO] No protein discovery files were found.\n'
                '(debug) ->\n'
                f'\tinput_path: {input_path}'
            )
            _logger.info(message)
            print(message)
            return

        #######################################################################
        # Load all discovered proteins into the default database.
        # 1. Ensure the desired database and tables exist.
        ORFaqsProteinQueryApi.create_workspace(workspace)
        ORFaqsProteinQueryApi._create_orfaqs_discovered_proteins_table()

        # 2. Connect the the workspace database.
        database_connection = PostgresDatabaseUtils.connect(
            _database_connection_options
        )

        # 3. Load all the result files into memory and write them to the
        # database table.
        for file_path in discovered_proteins_files:
            proteins_dataframe = PandasUtils.read_file_as_dataframe(file_path)
            ORFaqsProteinQueryApi._prepare_dataframe_for_table(
                proteins_dataframe
            )
            PostgresDatabaseUtils.insert_from_dataframe(
                connection=database_connection,
                table=_ORFaqsDiscoveredProteinsTableFactory.TABLE_NAME,
                dataframe=proteins_dataframe,
                insert_method='append',
            )

    @staticmethod
    def load_reference_proteins(
        workspace: str,
        input_path: (str | os.PathLike),
    ):
        """
        Reads the contents of the file or directory path provided and loads the
        data into the defined database. If the provided path is a directory, it
        must contain the file ORFaqsProteinDiscoveryUtils.exported_file_name()

        Attempts to load data are only made if data checks are passed.
        Otherwise, the process is aborted and no data is available.

        ---------
        Arguments
        ---------
        workspace (str):
            The name of the workspace session to load the discovered proteins.
            Workspaces persist until they are removed by the calling
            application.

        input_path (str | os.PathLike):
            A file or directory path to discovered proteins.
        """
        #######################################################################
        # Gather all discovery files.
        input_file_paths: list[os.PathLike] = []
        if DirectoryUtils.is_file(input_path):
            input_file_paths.append(input_path)
        elif DirectoryUtils.is_directory(input_path):
            # Grab all files from the directory.
            input_file_paths = DirectoryUtils.glob_files(
                input_path, recursive=True
            )

        fasta_files: list[os.PathLike] = []
        for file_path in input_file_paths:
            if FASTAUtils.is_fasta_file(file_path):
                fasta_files.append(file_path)

        if len(fasta_files) == 0:
            message = (
                '[INFO] No protein FASTA files were found.\n'
                '(debug) ->\n'
                f'\tinput_path: {input_path}'
            )
            _logger.info(message)
            print(message)
            return

        #######################################################################
        # Load all fasta proteins into the default database.
        # 1. Ensure the desired database and tables exist.
        ORFaqsProteinQueryApi.create_workspace(workspace)
        ORFaqsProteinQueryApi._create_orfaqs_reference_proteins_table()

        # 2. Connect the the workspace database.
        database_connection = PostgresDatabaseUtils.connect(
            _database_connection_options
        )

        # 3. Load all the result files into memory and write them to the
        # database table.
        for file_path in fasta_files:
            fasta_sequences = FASTAUtils.parse_file(file_path)
            protein_records: list[ORFaqsProteinRecord] = []
            for sequence in fasta_sequences:
                if sequence.sequence_type is not FASTASequenceType.AMINO_ACID:
                    continue

                protein_records.append(
                    ORFaqsProteinRecord(sequence.uid, sequence.sequence)
                )

            # Organize the results into a DataFrame object.
            records_dataframe = ORFaqsRecordUtils.orfaqs_records_to_dataframe(
                protein_records
            )
            PostgresDatabaseUtils.insert_from_dataframe(
                connection=database_connection,
                table=_ORFaqsReferenceProteinsTableFactory.TABLE_NAME,
                dataframe=records_dataframe,
                insert_method='append',
            )

    @staticmethod
    def export_proteins(
        workspace: str,
        file_path: (str | os.PathLike),
        export_format: _ExportFormatOptions = None,
        query_condition: str = None,
    ):
        if export_format is None:
            export_format = ORFaqsProteinQueryApi.default_export_format()

        ORFaqsProteinQueryApi.set_workspace(workspace)
        database_connection = PostgresDatabaseUtils.connect(
            _database_connection_options
        )

        query = f'SELECT * FROM {ORFaqsProteinQueryApi.table()}'
        if isinstance(query_condition, str):
            query += f' {query_condition}'

        results_dataframe = pd.read_sql(
            sql=query,
            con=database_connection,
        )
        ORFaqsProteinQueryApi._prepare_dataframe_for_export(results_dataframe)
        PandasUtils.export_dataframe(
            file_path=file_path,
            dataframe=results_dataframe,
            export_format=export_format,
        )
