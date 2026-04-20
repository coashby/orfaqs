"""
ORFaqs Protein Query common app classes, resources, and utility functions.
"""

import enum
import logging
import os
import pandas as pd
import pathlib
import typing

from orfaqs.modules.python.common.orfaqsapi import ORFaqsApi

from orfaqs.modules.python.proteindiscovery.orfaqsproteindiscovery import (
    ORFaqsProteinDiscoveryApi,
    RNAReadingFrame,
)
from orfaqs.modules.python.orfaqsrecords.orfaqsrecords import (
    ORFaqsDiscoveredProteinRecordKeys,
    ORFaqsProteinRecord,
    ORFaqsRecordUtils,
)

from orfaqs.lib.python.core.nucleotides import NucleotideUtils
from orfaqs.lib.python.utils.databaseutils import (
    BaseTable,
    PostgresDatabaseUtils,
    sqlalchemy,
    SqlAlchemyUtils,
)
from orfaqs.lib.python.utils.directoryutils import DirectoryUtils
from orfaqs.lib.python.utils.fastautils import (
    FASTASequenceType,
    FASTAUtils,
)
from orfaqs.lib.python.utils.pandasutils import (
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
        return f'The genomic sequence strand {strand_types_str} used during translation.'

    @staticmethod
    def reading_frame_comment() -> str:
        return 'The reading frame (1, 2, or 3) from which the protein was transcribed.'

    @staticmethod
    def rna_sequence_position_comment() -> str:
        return (
            'The position of the start codon in the RNA sequence '
            'from which the protein was transcribed.'
        )

    @staticmethod
    def protein_comment() -> str:
        return 'The protein represented by its single-letter amino acid abbreviations.'

    @staticmethod
    def protein_length_comment() -> str:
        return 'The number of amino acids in the protein sequence.'


class _ORFaqsProteinTableType(enum.Enum):
    DISCOVERED_PROTEINS = enum.auto()
    REFERENCE_PROTEINS = enum.auto()


class _ORFaqsDiscoveredProteinsTableFactory:
    """_ORFaqsDiscoveredProteinsTableUtil"""

    TABLE_NAME = 'discovered_proteins'

    @staticmethod
    def define_table(table_name: str = None) -> str:
        """
        Defines the discovered_proteins table when called iff the table
        definition does NOT exists in the base class.
        """
        if table_name is None:
            table_name = _ORFaqsDiscoveredProteinsTableFactory.TABLE_NAME
        if table_name in BaseTable.metadata.tables:
            return

        class ORFaqsDiscoveredProteinsTable(BaseTable):
            """ORFaqsDiscoveredProteinTable"""

            __tablename__ = table_name
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

        return table_name


class _ORFaqsReferenceProteinsTableFactory:
    """_ORFaqsReferenceProteinsTableFactory"""

    TABLE_NAME = 'reference_proteins'

    @staticmethod
    def define_table(table_name: str = None) -> str:
        """
        Defines the discovered_proteins table when called iff the table
        definition does NOT exists in the base class.
        """
        if table_name is None:
            table_name = _ORFaqsReferenceProteinsTableFactory.TABLE_NAME
        if table_name in BaseTable.metadata.tables:
            return

        class _ORFaqsReferenceProteinsTable(BaseTable):
            """_ORFaqsReferenceProteinsTable"""

            __tablename__ = table_name
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

        return table_name


class ORFaqsProteinQueryApi(ORFaqsApi):
    """ORFaqsProteinQueryApi"""

    _SESSION_DATA_WORKSPACES_KEY = 'workspaces'

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
    def default_table() -> str:
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
    def _create_orfaqs_discovered_proteins_table(
        table_name: str = None,
    ) -> str:
        table_name = _ORFaqsDiscoveredProteinsTableFactory.define_table(
            table_name
        )
        PostgresDatabaseUtils.create_table(
            BaseTable,
            _database_connection_options,
        )
        return table_name

    @staticmethod
    def _create_orfaqs_reference_proteins_table(table_name: str = None) -> str:
        table_name = _ORFaqsReferenceProteinsTableFactory.define_table(
            table_name
        )
        PostgresDatabaseUtils.create_table(
            BaseTable,
            _database_connection_options,
        )
        return table_name

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
    def _get_current_session_workspaces() -> dict:
        current_session_data: dict = {}
        session_data_object = ORFaqsProteinQueryApi._load_session_file()
        if session_data_object is not None:
            current_session_data: dict = session_data_object.session_data
        if (
            ORFaqsProteinQueryApi._SESSION_DATA_WORKSPACES_KEY
            not in current_session_data
        ):
            current_session_data[
                ORFaqsProteinQueryApi._SESSION_DATA_WORKSPACES_KEY
            ] = {}

        return current_session_data.get(
            ORFaqsProteinQueryApi._SESSION_DATA_WORKSPACES_KEY
        )

    @staticmethod
    def managed_workspaces() -> list[str]:
        session_workspaces = (
            ORFaqsProteinQueryApi._get_current_session_workspaces()
        )

        return sorted(list(session_workspaces.keys()))

    @staticmethod
    def managed_workspace_tables(workspace: str) -> list[str]:
        session_workspaces = (
            ORFaqsProteinQueryApi._get_current_session_workspaces()
        )
        workspace_tables: list[str] = []
        if workspace in session_workspaces:
            workspace_tables = session_workspaces[workspace]

        return sorted(workspace_tables)

    @staticmethod
    def _update_session_workspaces(workspaces_dict: dict[str, list]):
        new_session_data: dict = {
            ORFaqsProteinQueryApi._SESSION_DATA_WORKSPACES_KEY: workspaces_dict
        }
        ORFaqsProteinQueryApi._update_session_file(
            session_data=new_session_data,
        )

    @staticmethod
    def _add_record_to_session_data(
        workspace: str,
        table_name: str = None,
    ):
        workspaces_dict: dict[str, list] = (
            ORFaqsProteinQueryApi._get_current_session_workspaces()
        )
        if workspace not in workspaces_dict:
            workspaces_dict[workspace] = []

        if table_name is not None:
            workspace_tables: list = workspaces_dict[workspace]
            if table_name not in workspace_tables:
                workspace_tables.append(table_name)

            workspaces_dict[workspace] = sorted(workspace_tables)

        ORFaqsProteinQueryApi._update_session_workspaces(workspaces_dict)

    @staticmethod
    def _remove_record_from_session_data(
        workspace: str,
        table_name: str = None,
    ):
        workspaces_dict: dict[str, list] = (
            ORFaqsProteinQueryApi._get_current_session_workspaces()
        )
        if table_name is None:
            try:
                del workspaces_dict[workspace]
            except ValueError:
                return
        else:
            try:
                workspace_tables: list = workspaces_dict[workspace]
                workspace_tables.remove(table_name)
                workspaces_dict[workspace] = sorted(workspace_tables)
            except ValueError:
                return

        ORFaqsProteinQueryApi._update_session_workspaces(workspaces_dict)

    @staticmethod
    def create_workspace(workspace: str):
        ORFaqsProteinQueryApi.set_workspace(workspace)
        PostgresDatabaseUtils.create_database(_database_connection_options)
        ORFaqsProteinQueryApi._add_record_to_session_data(workspace=workspace)

    @staticmethod
    def _create_orfaqs_protein_table(
        workspace: str,
        table_name: str,
        table_type: _ORFaqsProteinTableType,
    ) -> str:
        ORFaqsProteinQueryApi.set_workspace(workspace)
        PostgresDatabaseUtils.create_database(_database_connection_options)
        if _ORFaqsProteinTableType.DISCOVERED_PROTEINS == table_type:
            table_name = (
                ORFaqsProteinQueryApi._create_orfaqs_discovered_proteins_table(
                    table_name
                )
            )
        elif _ORFaqsProteinTableType.REFERENCE_PROTEINS == table_type:
            table_name = (
                ORFaqsProteinQueryApi._create_orfaqs_reference_proteins_table(
                    table_name
                )
            )

        ORFaqsProteinQueryApi._add_record_to_session_data(
            workspace=workspace,
            table_name=table_name,
        )

        return table_name

    @staticmethod
    def remove_workspace(workspace: str):
        # Only remove managed workspaces.
        if workspace not in ORFaqsProteinQueryApi.managed_workspaces():
            return

        PostgresDatabaseUtils.drop_database(workspace)
        ORFaqsProteinQueryApi._remove_record_from_session_data(
            workspace=workspace
        )

    @staticmethod
    def remove_table(
        workspace: str,
        table_name: str = None,
    ):
        if table_name is None:
            table_name = ORFaqsProteinQueryApi.default_table()

        # Only remove managed workspace tables.
        managed_tables = ORFaqsProteinQueryApi.managed_workspace_tables(
            workspace
        )
        if table_name not in managed_tables:
            return
        ORFaqsProteinQueryApi.set_workspace(workspace)
        PostgresDatabaseUtils.drop_table(
            table=table_name,
            connection_options=_database_connection_options,
        )
        ORFaqsProteinQueryApi._remove_record_from_session_data(
            workspace=workspace,
            table_name=table_name,
        )

    @staticmethod
    def _load_discovered_proteins(
        workspace: str,
        table_name: str,
        discovered_proteins_files: list[pathlib.Path],
    ):
        #######################################################################
        # Load all discovered proteins into the default database.
        # 1. Ensure the desired database and tables exist.
        table_name = ORFaqsProteinQueryApi._create_orfaqs_protein_table(
            workspace=workspace,
            table_name=table_name,
            table_type=_ORFaqsProteinTableType.DISCOVERED_PROTEINS,
        )

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
                table=table_name,
                dataframe=proteins_dataframe,
                insert_method='append',
            )

    @staticmethod
    def _load_fasta_proteins(
        workspace: str,
        table_name: str,
        fasta_protein_files: list[pathlib.Path],
    ):
        #######################################################################
        # Load all fasta proteins into the default database.
        # 1. Ensure the desired database and tables exist.
        table_name = ORFaqsProteinQueryApi._create_orfaqs_protein_table(
            workspace=workspace,
            table_name=table_name,
            table_type=_ORFaqsProteinTableType.REFERENCE_PROTEINS,
        )

        # 2. Connect the the workspace database.
        database_connection = PostgresDatabaseUtils.connect(
            _database_connection_options
        )

        # 3. Load all the result files into memory and write them to the
        # database table.
        for file_path in fasta_protein_files:
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
                table=table_name,
                dataframe=records_dataframe,
                insert_method='append',
            )

    @staticmethod
    def load_proteins(
        workspace: str,
        table_name: str,
        input_path: (str | os.PathLike),
    ):
        """
        Reads the contents of the file or directory path provided and loads the
        data into the defined database.

        Attempts to load data are only made if data checks are passed.
        Otherwise, the process is aborted and no data is available.

        ---------
        Arguments
        ---------
        workspace (str):
            The name of the workspace session to load the discovered proteins.
            Workspaces persist until they are removed by the calling
            application.

        table (str):
            The name of the workspace session to load the discovered proteins.
            Workspaces persist until they are removed by the calling
            application.

        input_path (str | os.PathLike):
            A file or directory path to discovered proteins or FASTA files.
        """
        discovered_proteins_files = (
            ORFaqsProteinDiscoveryApi.find_discovered_protein_files(input_path)
        )
        if len(discovered_proteins_files) > 0:
            ORFaqsProteinQueryApi._load_discovered_proteins(
                workspace=workspace,
                table_name=table_name,
                discovered_proteins_files=discovered_proteins_files,
            )
            return

        fasta_protein_files = FASTAUtils.find_fasta_files(input_path)
        if len(fasta_protein_files) > 0:
            ORFaqsProteinQueryApi._load_fasta_proteins(
                workspace=workspace,
                table_name=table_name,
                fasta_protein_files=fasta_protein_files,
            )
            return

        message = (
            '[INFO] No protein files were found in the given '
            'file/directory path.\n'
            '(debug) ->\n'
            f'\tinput_path: {input_path}'
        )
        _logger.info(message)
        print(message)

    @staticmethod
    def export_proteins(
        workspace: str,
        table_name: str = None,
        export_path: (str | os.PathLike) = None,
        export_format: _ExportFormatOptions = None,
        query_condition: str = None,
    ):
        if table_name is None:
            table_name = ORFaqsProteinQueryApi.default_table()
        if export_path is None:
            export_path = './'
        if export_format is None:
            export_format = ORFaqsProteinQueryApi.default_export_format()

        ORFaqsProteinQueryApi.set_workspace(workspace)
        database_connection = PostgresDatabaseUtils.connect(
            _database_connection_options
        )

        query = f'SELECT * FROM {table_name}'
        if isinstance(query_condition, str):
            query += f' {query_condition}'

        results_dataframe = pd.read_sql(
            sql=query,
            con=database_connection,
        )
        ORFaqsProteinQueryApi._prepare_dataframe_for_export(results_dataframe)
        PandasUtils.export_dataframe(
            file_path=export_path,
            dataframe=results_dataframe,
            export_format=export_format,
        )
