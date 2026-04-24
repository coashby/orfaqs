#!/usr/bin/env python
"""
ORFaqs Protein Query
"""

import enum
import logging
import os
import typer

from typing import (
    Annotated,
    Callable,
)

from orfaqs.modules.python.common.orfaqscli import ORFaqsCli
from orfaqs.modules.python.orfaqsproteinquery.orfaqsproteinquery import (
    _ExportFormatOptions,
    ORFaqsProteinQueryApi,
)
from orfaqs.lib.python.utils.jsonutils import JsonUtils
from orfaqs.lib.python.utils.pandasutils import PandasUtils

_logger = logging.getLogger(__name__)
app = ORFaqsCli._new_app_typer()


class ORFaqsProteinQueryCli(ORFaqsCli):
    """ORFaqsProteinQuery"""

    _RESULTS_FILE_NAME = 'protein-query'
    _DEFAULT_WORKSPACE = ORFaqsProteinQueryApi.default_workspace()
    _DEFAULT_TABLE_NAME = ORFaqsProteinQueryApi.default_table()

    class Commands(enum.Enum):
        LOAD_PROTEINS = 'load-proteins'
        QUERY = 'query'
        EXPORT = 'export'
        REMOVE_TABLE = 'remove-table'
        SHOW = 'show'
        REMOVE_WORKSPACE = 'remove-workspace'

    @staticmethod
    def program_name() -> str:
        return 'ORFaqs Protein Query'

    @staticmethod
    def default_workspace() -> str:
        return ORFaqsProteinQueryCli._DEFAULT_WORKSPACE

    @staticmethod
    def default_table_name() -> str:
        return ORFaqsProteinQueryCli._DEFAULT_TABLE_NAME

    @staticmethod
    def _workspace_annotation(parameter_type: Callable) -> tuple:
        return (
            str,
            parameter_type(
                help=(
                    'The workspace (database) to use for the command. '
                    'Workspace names must be alphanumeric [Aa-Zz, 0-9]. The '
                    'underscore character "_" is permitted. '
                ),
            ),
        )

    @staticmethod
    def _table_name_annotation(parameter_type: Callable) -> tuple:
        return (
            str,
            parameter_type(
                help=('The table to use for the command.'),
            ),
        )

    @staticmethod
    def _query_annotation(parameter_type: Callable) -> tuple:
        return (
            str,
            parameter_type(
                help=(
                    'Query string for inspecting the protein sequence or set '
                    'of sequences specified in the input.'
                )
            ),
        )

    @staticmethod
    def _load_proteins(
        proteins: (str | os.PathLike),
        workspace: str = None,
        table: str = None,
        **kwargs,
    ):
        if workspace is None:
            workspace = ORFaqsProteinQueryCli._DEFAULT_WORKSPACE

        ORFaqsProteinQueryApi.load_proteins(
            workspace=workspace,
            table=table,
            input_path=proteins,
        )

    @app.command(Commands.LOAD_PROTEINS.value)
    @staticmethod
    def load_proteins(
        proteins: Annotated[
            str,
            typer.Argument(
                help=(
                    'A protein file, or directory containing protein data. '
                    'Data is loaded into the specified workspace and database '
                    'table.'
                )
            ),
        ],
        workspace: Annotated[_workspace_annotation(typer.Option)] = None,
        table: Annotated[_table_name_annotation(typer.Option)] = None,
        launch_json: Annotated[ORFaqsCli._launch_json_annotation()] = None,
        ctx: typer.Context = None,
    ):
        ui_kwargs = ORFaqsCli._process_ui_kwargs(**ctx.params)
        ORFaqsProteinQueryCli._load_proteins(**ui_kwargs)

    @staticmethod
    def _query_proteins(
        query: str = None,
        workspace: str = None,
        table: str = None,
        output_directory: str | os.PathLike = None,
        export_path: str | os.PathLike = None,
        job_id: str = None,
        export_format: _ExportFormatOptions = None,
        **kwargs,
    ):
        if workspace is None:
            workspace = ORFaqsProteinQueryCli.default_workspace()

        if export_path is None and export_format is None:
            query_result = ORFaqsProteinQueryApi.query_proteins(
                query,
                workspace=workspace,
                table=table,
            )
            PandasUtils.print_dataframe(query_result)
        else:
            if export_path is None:
                export_path = ORFaqsProteinQueryCli._create_output_directory(
                    output_directory=output_directory,
                    job_id=job_id,
                )

            ORFaqsProteinQueryApi.export_proteins(
                workspace=workspace,
                table=table,
                export_path=export_path,
                export_format=export_format,
                query=query,
            )

    @app.command(Commands.QUERY.value)
    @staticmethod
    def query_proteins(
        query: Annotated[_query_annotation(typer.Argument)] = None,
        workspace: Annotated[_workspace_annotation(typer.Option)] = None,
        table: Annotated[_table_name_annotation(typer.Option)] = None,
        output_directory: Annotated[
            ORFaqsCli._output_directory_annotation()
        ] = None,
        job_id: Annotated[ORFaqsCli._job_id_annotation()] = None,
        export_path: Annotated[ORFaqsCli._export_path_annotation()] = None,
        export_format: Annotated[
            ORFaqsCli._export_format_annotation(_ExportFormatOptions)
        ] = None,
        launch_json: Annotated[ORFaqsCli._launch_json_annotation()] = None,
        ctx: typer.Context = None,
    ):
        ui_kwargs = ORFaqsCli._process_ui_kwargs(**ctx.params)
        ORFaqsProteinQueryCli._query_proteins(**ui_kwargs)

    @staticmethod
    def _export_proteins(
        query: str = None,
        workspace: str = None,
        table: str = None,
        output_directory: str | os.PathLike = None,
        export_path: str | os.PathLike = None,
        job_id: str = None,
        export_format: _ExportFormatOptions = None,
        **kwargs,
    ):
        #######################################################################
        # Create the local output directory path
        output_directory_path = ORFaqsProteinQueryCli._create_output_directory(
            output_directory,
            job_id,
        )
        #######################################################################
        # Export proteins conditions...
        if export_path is None:
            export_path = (
                output_directory_path
                / ORFaqsProteinQueryApi.default_export_file_name()
            )

        ORFaqsProteinQueryApi.export_proteins(
            workspace=workspace,
            table=table,
            export_path=export_path,
            export_format=export_format,
            query=query,
        )

    @app.command(Commands.EXPORT.value)
    @staticmethod
    def export_proteins(
        workspace: Annotated[_workspace_annotation(typer.Option)] = None,
        table: Annotated[_table_name_annotation(typer.Option)] = None,
        output_directory: Annotated[
            ORFaqsCli._output_directory_annotation()
        ] = None,
        job_id: Annotated[ORFaqsCli._job_id_annotation()] = None,
        export_path: Annotated[ORFaqsCli._export_path_annotation()] = None,
        export_format: Annotated[
            ORFaqsCli._export_format_annotation(_ExportFormatOptions)
        ] = ORFaqsProteinQueryApi.default_export_format(),
        query: Annotated[_query_annotation(typer.Option)] = None,
        launch_json: Annotated[ORFaqsCli._launch_json_annotation()] = None,
        ctx: typer.Context = None,
    ):
        ui_kwargs = ORFaqsCli._process_ui_kwargs(**ctx.params)
        ORFaqsProteinQueryCli._export_proteins(**ui_kwargs)

    @staticmethod
    def _remove_workspace(
        workspace: str = None,
        all_workspaces: bool = False,
        **kwargs,
    ):
        if all_workspaces:
            ORFaqsProteinQueryApi.remove_all_managed_workspaces()
        else:
            ORFaqsProteinQueryApi.remove_workspace(workspace)

    @app.command(Commands.REMOVE_WORKSPACE.value)
    @staticmethod
    def remove_workspace(
        workspace: Annotated[_workspace_annotation(typer.Argument)] = None,
        all_workspaces: Annotated[
            bool,
            typer.Option(
                '--all',
                help=(
                    'Remove all managed workspaces in the database server.'
                    'specified workspace, or the default workspace if no workspace is specified.'
                ),
            ),
        ] = False,
        ctx: typer.Context = None,
    ):
        if 'all' in ctx.params:
            ctx.params['all_workspaces'] = ctx.pop('all', False)

        ui_kwargs = ORFaqsCli._process_ui_kwargs(**ctx.params)
        ORFaqsProteinQueryCli._remove_workspace(**ui_kwargs)

    @staticmethod
    def _remove_table(
        workspace: str,
        table: str = None,
        all_tables: bool = False,
        **kwargs,
    ):
        if all_tables:
            ORFaqsProteinQueryApi.remove_all_managed_tables(workspace)
        else:
            ORFaqsProteinQueryApi.remove_table(
                workspace=workspace,
                table=table,
            )

    @app.command(Commands.REMOVE_TABLE.value)
    @staticmethod
    def remove_table(
        workspace: Annotated[_workspace_annotation(typer.Argument)] = None,
        table: Annotated[_table_name_annotation(typer.Argument)] = None,
        all_tables: Annotated[
            bool,
            typer.Option(
                '--all',
                help=(
                    'Remove all managed tables in the specified workspace, or '
                    'the default workspace if no workspace is specified.'
                ),
            ),
        ] = False,
        ctx: typer.Context = None,
    ):
        if 'all' in ctx.params:
            ctx.params['all_tables'] = ctx.pop('all', False)

        ui_kwargs = ORFaqsCli._process_ui_kwargs(**ctx.params)
        ORFaqsProteinQueryCli._remove_table(**ui_kwargs)

    @staticmethod
    def _show(
        workspace: str = None,
        workspaces: bool = False,
        all_workspaces_and_tables: bool = False,
        **kwargs,
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

    @app.command(Commands.SHOW.value)
    @staticmethod
    def show(
        workspace: Annotated[
            str,
            typer.Option(
                help=(
                    'Display all managed tables in the specified workspace.'
                ),
            ),
        ] = None,
        workspaces: Annotated[
            bool,
            typer.Option(
                '--workspaces',
                help=('Display all managed workspaces.'),
            ),
        ] = False,
        all_workspaces_and_tables: Annotated[
            bool,
            typer.Option(
                '--all',
                help=('Display all managed workspaces and tables.'),
            ),
        ] = False,
        ctx: typer.Context = None,
    ):
        if 'all' in ctx.params:
            ctx.params['all_workspaces_and_tables'] = ctx.pop('all', False)

        ui_kwargs = ORFaqsCli._process_ui_kwargs(**ctx.params)
        ORFaqsProteinQueryCli._show(**ui_kwargs)

    @staticmethod
    def _run_cli_command(
        cli_command: Commands,
        ctx: typer.Context,
    ):
        ctx.params[ORFaqsCli._KWARGS_CTX_OBJECT_KEY] = ctx
        if ORFaqsProteinQueryCli.Commands.LOAD_PROTEINS == cli_command:
            ORFaqsProteinQueryCli.load_proteins(**ctx.params)
        elif ORFaqsProteinQueryCli.Commands.QUERY == cli_command:
            ORFaqsProteinQueryCli.query_proteins(**ctx.params)
        elif ORFaqsProteinQueryCli.Commands.EXPORT == cli_command:
            ORFaqsProteinQueryCli.export_proteins(**ctx.params)
        elif ORFaqsProteinQueryCli.Commands.REMOVE_WORKSPACE == cli_command:
            ORFaqsProteinQueryCli.remove_workspace(**ctx.params)
        elif ORFaqsProteinQueryCli.Commands.REMOVE_TABLE == cli_command:
            ORFaqsProteinQueryCli.remove_table(**ctx.params)
        elif ORFaqsProteinQueryCli.Commands.SHOW == cli_command:
            ORFaqsProteinQueryCli.show(**ctx.params)

    @app.callback(invoke_without_command=True)
    @staticmethod
    def cli(
        launch_json: Annotated[ORFaqsCli._launch_json_annotation()] = None,
        ctx: typer.Context = None,
    ):
        ui_kwargs = ORFaqsCli._process_ui_kwargs(**ctx.params)
        number_of_commands = len(ui_kwargs)
        number_commands_executed = 0
        for cli_command in ORFaqsProteinQueryCli.Commands:
            command_arg = ORFaqsCli._cli_params_to_args(cli_command)
            if command_arg in ui_kwargs:
                # Update the typer Context.
                ctx.params = ui_kwargs.get(command_arg, {})
                ORFaqsProteinQueryCli._run_cli_command(
                    cli_command=cli_command,
                    ctx=ctx,
                )
                number_commands_executed += 1
            if number_commands_executed == number_of_commands:
                break

    @staticmethod
    def _run():
        app()


if __name__ == '__main__':
    ORFaqsProteinQueryCli.run()
