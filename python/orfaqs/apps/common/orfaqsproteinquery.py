"""
ORFaqs Protein Query common app classes, resources, and utility functions.
"""

import logging
import os
import pandas as pd
import psycopg2
import typing

from orfaqs.apps.common.orfaqsproteindiscovery import (
    ORFaqsProteinDiscoveryUtils,
    ORFaqsProteinDiscoveryRecordKeys,
    RNAReadingFrame,
)

from orfaqs.lib.core.nucleotides import NucleotideUtils
from orfaqs.lib.utils.databaseutils import (
    BaseTable,
    PostgresDatabaseUtils,
    sqlalchemy,
    SqlAlchemyUtils,
)
from orfaqs.lib.utils.directoryutils import DirectoryUtils
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


class ORFaqsDiscoveredProteinTableSchema(ORFaqsProteinDiscoveryRecordKeys):
    UID_KEY = 'uid'

    @staticmethod
    def columns() -> list[str]:
        return [
            ORFaqsDiscoveredProteinTableSchema.UID_KEY,
            ORFaqsDiscoveredProteinTableSchema.ACCESSION_NUMBER_KEY,
            ORFaqsDiscoveredProteinTableSchema.STRAND_TYPE_KEY,
            ORFaqsDiscoveredProteinTableSchema.READING_FRAME_KEY,
            ORFaqsDiscoveredProteinTableSchema.RNA_SEQUENCE_POSITION_KEY,
            ORFaqsDiscoveredProteinTableSchema.PROTEIN_KEY,
            ORFaqsDiscoveredProteinTableSchema.PROTEIN_LENGTH_KEY,
        ]


class _ORFaqsDiscoveredProteinsTableFactory:
    """_ORFaqsDiscoveredProteinsTableUtil"""

    TABLE_NAME = 'discovered_proteins'

    @staticmethod
    def _define_table():
        """
        Defines the discovered_proteins table when called iff the table
        definition does NOT exists in the base class.
        """
        table_name = 'discovered_proteins'
        if table_name in BaseTable.metadata.tables:
            return

        class ORFaqsDiscoveredProteinsTable(BaseTable):
            """ORFaqsDiscoveredProteinTable"""

            __tablename__ = _ORFaqsDiscoveredProteinsTableFactory.TABLE_NAME
            ACCESSION_NUMBER_MAX_CHAR_LENGTH = 64
            STRAND_TYPE_MAX_CHAR_LENGTH = 32
            UID_MAX_CHAR_LENGTH = 64

            @staticmethod
            def _reading_frame_check_constraint() -> str:
                reading_frame_values = [
                    reading_frame.value for reading_frame in RNAReadingFrame
                ]
                return f'reading_frame IN {tuple(reading_frame_values)}'

            @staticmethod
            def uid_comment() -> str:
                return 'The unique id of the entry.'

            @staticmethod
            def accession_number_comment() -> str:
                return (
                    'The accession number assigned to the sequence from '
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

            uid = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.UID_KEY,
                sqlalchemy.String(UID_MAX_CHAR_LENGTH),
                primary_key=True,
                comment=uid_comment(),
            )

            accession_number = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.ACCESSION_NUMBER_KEY,
                sqlalchemy.String(ACCESSION_NUMBER_MAX_CHAR_LENGTH),
                comment=accession_number_comment(),
            )

            strand_type = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.STRAND_TYPE_KEY,
                sqlalchemy.String(STRAND_TYPE_MAX_CHAR_LENGTH),
                nullable=False,
                comment=strand_type_comment(),
            )

            reading_frame = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.READING_FRAME_KEY,
                sqlalchemy.Integer,
                nullable=False,
                comment=reading_frame_comment(),
            )

            rna_sequence_position = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.RNA_SEQUENCE_POSITION_KEY,
                sqlalchemy.Integer,
                sqlalchemy.CheckConstraint(_reading_frame_check_constraint()),
                nullable=False,
                comment=rna_sequence_position_comment(),
            )

            protein = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.PROTEIN_KEY,
                sqlalchemy.String,
                nullable=False,
                comment=protein_comment(),
            )

            protein_length = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinTableSchema.PROTEIN_LENGTH_KEY,
                sqlalchemy.Integer,
                nullable=False,
                comment=protein_length_comment(),
            )


class ORFaqsProteinQueryUtils:
    """ORFaqsProteinQueryUtils"""

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
    def default_database():
        return _ORFAQS_PROTEIN_QUERY_DATABASE

    @staticmethod
    def set_database(database: str):
        global _database_connection_options
        _database_connection_options.database = database

    @staticmethod
    def set_workspace(workspace: str):
        ORFaqsProteinQueryUtils.set_database(workspace)

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
        _ORFaqsDiscoveredProteinsTableFactory._define_table()
        PostgresDatabaseUtils.create_table(
            BaseTable,
            _database_connection_options,
        )

    @staticmethod
    def _create_uid(
        accession_number: str,
        strand_type: str,
        reading_frame: int,
        rna_sequence_position: int,
        protein_length: str,
    ) -> str:
        uid_str = ':'.join(
            [
                accession_number,
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
            accession_number = row[
                ORFaqsDiscoveredProteinTableSchema.ACCESSION_NUMBER_KEY
            ]
            strand_type = row[
                ORFaqsDiscoveredProteinTableSchema.STRAND_TYPE_KEY
            ]
            reading_frame = row[
                ORFaqsDiscoveredProteinTableSchema.READING_FRAME_KEY
            ]
            rna_sequence_position = row[
                ORFaqsDiscoveredProteinTableSchema.RNA_SEQUENCE_POSITION_KEY
            ]
            protein_length = row[
                ORFaqsDiscoveredProteinTableSchema.PROTEIN_LENGTH_KEY
            ]
            return ORFaqsProteinQueryUtils._create_uid(
                accession_number=accession_number,
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
        ORFaqsProteinQueryUtils._add_database_primary_key_column(
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
        ORFaqsProteinQueryUtils.set_workspace(workspace)
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
        discovery_files: list[os.PathLike] = []
        if DirectoryUtils.is_file(input_path):
            discovery_files.append(input_path)
        elif DirectoryUtils.is_directory(input_path):
            # Grab all files from the directory. In discoveries involving
            # multiple genes, discovery files will be organized in
            # subdirectories.
            file_paths = DirectoryUtils.glob_files(input_path, recursive=True)
            discovery_file_name = (
                ORFaqsProteinDiscoveryUtils.exported_file_name()
            )
            for file_path in file_paths:
                # Check that contents are indeed, formatted results from the
                # ORFaqsProteinDiscoveryUtils class.
                file_type = file_path.suffix
                expected_file_name_ending = f'{discovery_file_name}{file_type}'
                if (
                    (expected_file_name_ending in file_path.name)
                    and ORFaqsProteinDiscoveryUtils.is_valid_discovered_proteins_file(
                        file_path
                    )
                ):
                    discovery_files.append(file_path)

        if len(discovery_files) == 0:
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
        ORFaqsProteinQueryUtils.create_workspace(workspace)
        ORFaqsProteinQueryUtils._create_orfaqs_discovered_proteins_table()

        # 2. Connect the the workspace database.
        database_connection = PostgresDatabaseUtils.connect(
            _database_connection_options
        )

        # 3. Load all the result files into memory.
        for file_path in discovery_files:
            proteins_dataframe = PandasUtils.read_file_as_dataframe(file_path)
            ORFaqsProteinQueryUtils._prepare_dataframe_for_table(
                proteins_dataframe
            )
            try:
                proteins_dataframe.to_sql(
                    name=_ORFaqsDiscoveredProteinsTableFactory.TABLE_NAME,
                    con=database_connection,
                    if_exists='append',
                    index=False,
                    method=SqlAlchemyUtils.psql_insert_copy,
                )
            except psycopg2.errors.UniqueViolation:
                message = (
                    '[INFO] Duplicate data found in the current '
                    f'DataFrame object loaded from: {file_path}.'
                )
                _logger.info(message)
                print(message)

    @staticmethod
    def export_proteins(
        workspace: str,
        file_path: (str | os.PathLike),
        export_format: _ExportFormatOptions,
        query_condition: str = None,
    ):
        ORFaqsProteinQueryUtils.set_workspace(workspace)
        database_connection = PostgresDatabaseUtils.connect(
            _database_connection_options
        )

        query = f'SELECT * FROM {ORFaqsProteinQueryUtils.table()}'
        if isinstance(query_condition, str):
            query += f' {query_condition}'

        results_dataframe = pd.read_sql(
            sql=query,
            con=database_connection,
        )
        ORFaqsProteinQueryUtils._prepare_dataframe_for_export(
            results_dataframe
        )
        PandasUtils.export_dataframe(
            file_path=file_path,
            dataframe=results_dataframe,
            export_format=export_format,
        )
