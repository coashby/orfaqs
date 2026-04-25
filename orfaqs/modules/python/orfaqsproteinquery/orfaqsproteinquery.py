"""
ORFaqs Protein Query common app classes, resources, and utility functions.
"""

import logging
import os
import pandas as pd
import pathlib
import typing

from orfaqs.modules.python.common.orfaqsapi import ORFaqsApi

from orfaqs.modules.python.orfaqsproteindiscovery.orfaqsproteindiscovery import (
    ORFaqsProteinDiscoveryApi,
)
from orfaqs.modules.python.orfaqsproteinquery.databases.orfaqsproteintableutils import (
    ORFaqsProteinTableType,
)
from orfaqs.modules.python.orfaqsproteinquery.databases.models.orfaqsdiscoveredproteinstable import (
    ORFaqsDiscoveredProteinsTableFactory,
    ORFaqsDiscoveredProteinsTableSchema,
)
from orfaqs.modules.python.orfaqsproteinquery.databases.models.orfaqsreferenceproteinstable import (
    ORFaqsReferenceProteinsTableFactory,
)
from orfaqs.modules.python.orfaqsrecords.orfaqsrecords import (
    ORFaqsProteinRecord,
    ORFaqsRecordUtils,
)

from orfaqs.lib.python.utils.databaseutils import (
    BaseTable,
    PostgresDatabaseUtils,
)
from orfaqs.lib.python.utils.directoryutils import DirectoryUtils
from orfaqs.lib.python.utils.fastautils import (
    FASTASequenceType,
    FASTAUtils,
)
from orfaqs.lib.python.utils.jsonutils import JsonUtils
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

_ORFAQS_PROTEIN_QUERY_DEFAULT_DATABASE = 'orfaqs_protein_query'
_ORFAQS_PROTEIN_QUERY_DEFAULT_TABLE = 'proteins'
# Global database connection options enable user defined options to persist
# during the entirety of a session.
_database_connection_options = (
    PostgresDatabaseUtils.default_database_connection_options()
)
_database_connection_options.database = _ORFAQS_PROTEIN_QUERY_DEFAULT_DATABASE


class ORFaqsProteinQueryApi(ORFaqsApi):
    """ORFaqsProteinQueryApi"""

    _SESSION_DATA_WORKSPACES_KEY = 'workspaces'

    @staticmethod
    def _set_database(database: str):
        global _database_connection_options
        _database_connection_options.database = database

    @staticmethod
    def _configure_database(
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
        table: str = None,
    ) -> str:
        table = ORFaqsDiscoveredProteinsTableFactory.define_table(table)
        PostgresDatabaseUtils.create_table(
            BaseTable,
            _database_connection_options,
        )
        return table

    @staticmethod
    def _create_orfaqs_reference_proteins_table(table: str = None) -> str:
        table = ORFaqsReferenceProteinsTableFactory.define_table(table)
        PostgresDatabaseUtils.create_table(
            BaseTable,
            _database_connection_options,
        )
        return table

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
    def _add_discovered_proteins_table_primary_key_column(
        proteins_dataframe: pd.DataFrame,
    ):
        def create_uid(row: pd.Series):
            uid = row[ORFaqsDiscoveredProteinsTableSchema.SOURCE_UID_KEY]
            strand_type = row[
                ORFaqsDiscoveredProteinsTableSchema.STRAND_TYPE_KEY
            ]
            reading_frame = row[
                ORFaqsDiscoveredProteinsTableSchema.READING_FRAME_KEY
            ]
            rna_sequence_position = row[
                ORFaqsDiscoveredProteinsTableSchema.DNA_SEQUENCE_POSITION_KEY
            ]
            protein_length = row[
                ORFaqsDiscoveredProteinsTableSchema.PROTEIN_LENGTH_KEY
            ]
            return ORFaqsProteinQueryApi._create_uid(
                uid=uid,
                strand_type=strand_type,
                reading_frame=reading_frame,
                rna_sequence_position=rna_sequence_position,
                protein_length=protein_length,
            )

        proteins_dataframe[ORFaqsDiscoveredProteinsTableSchema.UID_KEY] = (
            proteins_dataframe.apply(create_uid, axis='columns')
        )

    @staticmethod
    def _prepare_dataframe_for_discovered_proteins_table(
        proteins_dataframe: pd.DataFrame,
    ):
        # Create a primary key
        ORFaqsProteinQueryApi._add_discovered_proteins_table_primary_key_column(
            proteins_dataframe
        )
        expected_columns = ORFaqsDiscoveredProteinsTableSchema.columns()
        drop_columns = proteins_dataframe.columns.difference(expected_columns)
        proteins_dataframe.drop(drop_columns, axis='columns', inplace=True)

    @staticmethod
    def _prepare_dataframe_for_export(proteins_dataframe: pd.DataFrame):
        proteins_dataframe.drop(
            columns=[ORFaqsDiscoveredProteinsTableSchema.UID_KEY],
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
        table: str = None,
    ):
        workspaces_dict: dict[str, list] = (
            ORFaqsProteinQueryApi._get_current_session_workspaces()
        )
        if workspace not in workspaces_dict:
            workspaces_dict[workspace] = []

        if table is not None:
            workspace_tables: list = workspaces_dict[workspace]
            if table not in workspace_tables:
                workspace_tables.append(table)

            workspaces_dict[workspace] = sorted(workspace_tables)

        ORFaqsProteinQueryApi._update_session_workspaces(workspaces_dict)

    @staticmethod
    def _remove_record_from_session_data(
        workspace: str,
        table: str = None,
    ):
        workspaces_dict: dict[str, list] = (
            ORFaqsProteinQueryApi._get_current_session_workspaces()
        )
        if table is None:
            try:
                del workspaces_dict[workspace]
            except ValueError:
                return
        else:
            try:
                workspace_tables: list = workspaces_dict[workspace]
                workspace_tables.remove(table)
                workspaces_dict[workspace] = sorted(workspace_tables)
            except ValueError:
                return

        ORFaqsProteinQueryApi._update_session_workspaces(workspaces_dict)

    @staticmethod
    def _create_orfaqs_protein_table(
        workspace: str,
        table: str,
        table_type: ORFaqsProteinTableType,
    ) -> str:
        ORFaqsProteinQueryApi.set_workspace(workspace)
        PostgresDatabaseUtils.create_database(_database_connection_options)
        if ORFaqsProteinTableType.DISCOVERED_PROTEINS == table_type:
            table = (
                ORFaqsProteinQueryApi._create_orfaqs_discovered_proteins_table(
                    table
                )
            )
        elif ORFaqsProteinTableType.REFERENCE_PROTEINS == table_type:
            table = (
                ORFaqsProteinQueryApi._create_orfaqs_reference_proteins_table(
                    table
                )
            )

        ORFaqsProteinQueryApi._add_record_to_session_data(
            workspace=workspace,
            table=table,
        )

        return table

    @staticmethod
    def _load_discovered_proteins(
        workspace: str,
        table: str,
        discovered_proteins_files: list[pathlib.Path],
    ):
        #######################################################################
        # Load all discovered proteins into the default database.
        # 1. Ensure the desired database and tables exist.
        table = ORFaqsProteinQueryApi._create_orfaqs_protein_table(
            workspace=workspace,
            table=table,
            table_type=ORFaqsProteinTableType.DISCOVERED_PROTEINS,
        )

        # 2. Connect the the workspace database.
        database_connection = PostgresDatabaseUtils.connect(
            _database_connection_options
        )

        # 3. Load all the result files into memory and write them to the
        # database table.
        for file_path in discovered_proteins_files:
            proteins_dataframe = PandasUtils.read_file_as_dataframe(file_path)
            ORFaqsProteinQueryApi._prepare_dataframe_for_discovered_proteins_table(
                proteins_dataframe
            )
            PostgresDatabaseUtils.insert_from_dataframe(
                connection=database_connection,
                table=table,
                dataframe=proteins_dataframe,
                insert_method='append',
            )

    @staticmethod
    def _load_fasta_proteins(
        workspace: str,
        table: str,
        fasta_protein_files: list[pathlib.Path],
    ):
        #######################################################################
        # Load all fasta proteins into the default database.
        # 1. Ensure the desired database and tables exist.
        table = ORFaqsProteinQueryApi._create_orfaqs_protein_table(
            workspace=workspace,
            table=table,
            table_type=ORFaqsProteinTableType.REFERENCE_PROTEINS,
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
                table=table,
                dataframe=records_dataframe,
                insert_method='append',
            )

    @staticmethod
    def _query_proteins(
        query: str = None,
        workspace: str = None,
        table: str = None,
    ) -> pd.DataFrame:
        if query is None and table is None:
            message = (
                '[ERROR] A valid query OR a table name must be specified.'
            )
            _logger.error(message)
            raise ValueError(message)

        if workspace is None:
            workspace = ORFaqsProteinQueryApi.default_workspace()
        ORFaqsProteinQueryApi.set_workspace(workspace)
        database_connection = PostgresDatabaseUtils.connect(
            _database_connection_options
        )

        if query is None:
            query = ''

        if 'select' not in query.lower() and table is not None:
            query = f'SELECT * FROM {table} {query}'

        query = query.rstrip()
        results_dataframe = pd.read_sql(
            sql=query,
            con=database_connection,
        )

        return results_dataframe

    @staticmethod
    def _export_dataframe(
        file_path: (str | os.PathLike),
        dataframe: pd.DataFrame,
        export_format: _ExportFormatOptions = None,
    ):
        ORFaqsProteinQueryApi._prepare_dataframe_for_export(dataframe)
        PandasUtils.export_dataframe(
            file_path=file_path,
            dataframe=dataframe,
            export_format=export_format,
        )

    @staticmethod
    def _export_proteins(
        workspace: str,
        table: str = None,
        export_path: (str | os.PathLike) = None,
        export_format: _ExportFormatOptions = None,
        query: str = None,
    ):
        if export_path is None:
            export_path = './'
        if export_format is None:
            export_format = ORFaqsProteinQueryApi.default_export_format()

        results_dataframe = ORFaqsProteinQueryApi._query_proteins(
            query,
            workspace=workspace,
            table=table,
        )

        export_file_path = export_path
        # If an explicit file path is not specified, use the table name OR
        # query as the file name for the export.
        if not DirectoryUtils.has_suffix(export_file_path):
            if workspace is not None:
                export_file_path = export_file_path / workspace

            if table is not None:
                export_file_path = export_file_path / table

            file_name = table
            if query is not None:
                file_name = DirectoryUtils.sanitize_path(
                    query.replace(' ', '_')
                )

            export_file_path = export_file_path / file_name

        ORFaqsProteinQueryApi._export_dataframe(
            file_path=export_file_path,
            dataframe=results_dataframe,
            export_format=export_format,
        )

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
        return _ORFAQS_PROTEIN_QUERY_DEFAULT_DATABASE

    @staticmethod
    def default_workspace():
        return ORFaqsProteinQueryApi.default_database()

    @staticmethod
    def set_workspace(workspace: str):
        ORFaqsProteinQueryApi._set_database(workspace)

    @staticmethod
    def default_table() -> str:
        return _ORFAQS_PROTEIN_QUERY_DEFAULT_TABLE

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
    def all_managed_workspaces_and_tables() -> dict[str, list[str]]:
        workspaces_tables_map: dict[str, list[str]] = {}
        for workspace in ORFaqsProteinQueryApi.managed_workspaces():
            workspaces_tables_map[workspace] = (
                ORFaqsProteinQueryApi.managed_workspace_tables(workspace)
            )

        return workspaces_tables_map

    @staticmethod
    def show_managed_workspaces(
        workspace: str = None,
        workspaces: bool = False,
        all_workspaces_and_tables: bool = False,
    ):
        results: dict[str, list[str]] = {}
        if all_workspaces_and_tables:
            results = ORFaqsProteinQueryApi.all_managed_workspaces_and_tables()
        elif workspaces:
            results['workspaces'] = ORFaqsProteinQueryApi.managed_workspaces()
        elif workspace is not None:
            results[workspace] = (
                ORFaqsProteinQueryApi.managed_workspace_tables(workspace)
            )

        if results is not None:
            results_json = JsonUtils.as_json_string(results)
            print(results_json)

    @staticmethod
    def create_workspace(workspace: str):
        ORFaqsProteinQueryApi.set_workspace(workspace)
        PostgresDatabaseUtils.create_database(_database_connection_options)
        ORFaqsProteinQueryApi._add_record_to_session_data(workspace=workspace)

    @staticmethod
    def _remove_workspace(workspace: str):
        # Only remove managed workspaces.
        if workspace not in ORFaqsProteinQueryApi.managed_workspaces():
            return

        PostgresDatabaseUtils.drop_database(workspace)
        ORFaqsProteinQueryApi._remove_record_from_session_data(
            workspace=workspace
        )

    @staticmethod
    def _remove_all_managed_workspaces():
        for workspace in ORFaqsProteinQueryApi.managed_workspaces():
            ORFaqsProteinQueryApi._remove_workspace(workspace)

    @staticmethod
    def remove_workspaces(
        workspaces: str | list[str] = None,
        all_workspaces: bool = False,
    ):
        if all_workspaces:
            ORFaqsProteinQueryApi._remove_all_managed_workspaces()
        else:
            if isinstance(workspaces, str):
                workspaces = workspaces.split(' ')

            for workspace in workspaces:
                ORFaqsProteinQueryApi._remove_workspace(workspace)

    @staticmethod
    def _remove_table(
        workspace: str,
        table: str = None,
    ):
        if table is None:
            table = ORFaqsProteinQueryApi.default_table()

        # Only remove managed workspace tables.
        managed_tables = ORFaqsProteinQueryApi.managed_workspace_tables(
            workspace
        )
        if table not in managed_tables:
            return
        ORFaqsProteinQueryApi.set_workspace(workspace)
        PostgresDatabaseUtils.drop_table(
            table=table,
            connection_options=_database_connection_options,
        )
        ORFaqsProteinQueryApi._remove_record_from_session_data(
            workspace=workspace,
            table=table,
        )

    @staticmethod
    def _remove_all_managed_tables(workspace: str):
        for table in ORFaqsProteinQueryApi.managed_workspace_tables(workspace):
            ORFaqsProteinQueryApi._remove_table(
                workspace=workspace,
                table=table,
            )

    @staticmethod
    def remove_tables(
        workspace: str,
        tables: str | list[str] = None,
        all_tables: bool = False,
    ):
        if all_tables:
            ORFaqsProteinQueryApi._remove_all_managed_tables(workspace)
        else:
            if isinstance(tables, str):
                tables = tables.split(' ')

            for table in tables:
                ORFaqsProteinQueryApi._remove_table(
                    workspace=workspace,
                    table=table,
                )

    @staticmethod
    def load_proteins(
        proteins: (str | os.PathLike),
        workspace: str = None,
        table: str = None,
        **kwargs,
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

        proteins (str | os.PathLike):
            A file or directory path to discovered proteins or FASTA files.
        """
        if workspace is None:
            workspace = ORFaqsProteinQueryApi.default_workspace()
        if table is None:
            table = ORFaqsProteinQueryApi.default_table()
        discovered_proteins_files = (
            ORFaqsProteinDiscoveryApi.find_discovered_protein_files(proteins)
        )
        if len(discovered_proteins_files) > 0:
            ORFaqsProteinQueryApi._load_discovered_proteins(
                workspace=workspace,
                table=table,
                discovered_proteins_files=discovered_proteins_files,
            )
            return

        fasta_protein_files = FASTAUtils.find_fasta_files(proteins)
        if len(fasta_protein_files) > 0:
            ORFaqsProteinQueryApi._load_fasta_proteins(
                workspace=workspace,
                table=table,
                fasta_protein_files=fasta_protein_files,
            )
            return

        message = (
            '[INFO] No protein files were found in the given '
            'file/directory path.\n'
            '(debug) ->\n'
            f'\tinput_path: {proteins}'
        )
        _logger.info(message)
        print(message)

    @staticmethod
    def query_proteins(
        query: str = None,
        workspace: str = None,
        table: str = None,
        export_path: str | os.PathLike = None,
        export_format: _ExportFormatOptions = None,
    ) -> pd.DataFrame | None:
        query_result: pd.DataFrame = None
        if export_path is None and export_format is None:
            query_result = ORFaqsProteinQueryApi._query_proteins(
                query,
                workspace=workspace,
                table=table,
            )
        else:
            ORFaqsProteinQueryApi.export_proteins(
                workspace=workspace,
                table=table,
                export_path=export_path,
                export_format=export_format,
                query=query,
            )

        return query_result

    @staticmethod
    def export_proteins(
        workspace: str = None,
        table: str = None,
        export_path: (str | os.PathLike) = None,
        export_format: _ExportFormatOptions = None,
        query: str = None,
    ):
        if workspace is None:
            workspace = ORFaqsProteinQueryApi.default_workspace()

        export_tables: list[str] = [table]
        if table is None and 'from' not in query.lower():
            # Export all data from all tables within the workspace.
            export_tables = ORFaqsProteinQueryApi.managed_workspace_tables(
                workspace
            )

        for table in export_tables:
            ORFaqsProteinQueryApi._export_proteins(
                workspace=workspace,
                table=table,
                export_path=export_path,
                export_format=export_format,
                query=query,
            )
