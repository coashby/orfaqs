#!/usr/bin/env python

import logging
import multiprocessing
import os
import typer

from typing import Annotated

from orfaqs.modules.python.common.orfaqscli import (
    ORFaqsCli,
    ORFaqsCliExitCodes,
)
from orfaqs.modules.python.orfaqsproteindiscovery.orfaqsproteindiscovery import (
    _ExportFormatOptions,
    ORFaqsProteinDiscoveryApi,
)

from orfaqs.lib.python.utils.directoryutils import DirectoryUtils


_logger = logging.getLogger(__name__)
app = ORFaqsCli._new_app_typer()


class ORFaqsProteinDiscoveryCli(ORFaqsCli):
    """ORFaqsProteinDiscovery"""

    @staticmethod
    def program_name():
        return 'ORFaqs Protein Discovery'

    @staticmethod
    def _discover_proteins(
        input_sequence: str | os.PathLike,
        uid: str = None,
        include_reverse_complement: bool = True,
        output_directory: str | os.PathLike = None,
        job_id: str = None,
        export_format: str = None,
        enable_gpu: bool = False,
    ):
        #######################################################################
        # Create the local output directory path
        if output_directory is None:
            output_directory = './'

        output_directory = DirectoryUtils.make_path_object(output_directory)
        output_directory = output_directory.joinpath(
            ORFaqsProteinDiscoveryCli.default_output_directory()
        )
        output_directory = DirectoryUtils.make_path_object(output_directory)
        if isinstance(job_id, str):
            output_directory = output_directory.joinpath(job_id)

        DirectoryUtils.mkdir_path(output_directory)

        if DirectoryUtils.is_file(input_sequence):
            # Try processing as a FASTA file
            ORFaqsProteinDiscoveryApi.process_fasta_file(
                fasta_file_path=input_sequence,
                include_reverse_complement=include_reverse_complement,
                export_format=export_format,
                output_directory=output_directory,
                use_gpu=enable_gpu,
            )
        else:
            # Try processing as an sequence string.
            ORFaqsProteinDiscoveryApi.process_genomic_sequence(
                genomic_sequence=input_sequence,
                uid=uid,
                include_reverse_complement=include_reverse_complement,
                export_format=export_format,
                output_directory=output_directory,
                use_gpu=enable_gpu,
            )

    @app.callback(invoke_without_command=True)
    @staticmethod
    def cli(
        input_sequence: Annotated[
            str,
            typer.Argument(
                help=(
                    'A sequence file or genomic sequence string. Currently, '
                    'only FASTA files are supported.'
                )
            ),
        ] = None,
        uid: Annotated[
            str,
            typer.Option(
                help=(
                    'Defines the accession number associated with the input'
                    'sequence string. If the input sequence is a file, this '
                    'option is ignored.'
                )
            ),
        ] = None,
        output_directory: Annotated[
            ORFaqsCli._output_directory_annotation()
        ] = None,
        job_id: Annotated[ORFaqsCli._job_id_annotation()] = None,
        export_format: Annotated[
            ORFaqsCli._export_format_annotation(_ExportFormatOptions)
        ] = ORFaqsProteinDiscoveryApi.default_export_format(),
        include_reverse_complement: Annotated[
            bool,
            typer.Option(
                '--include-reverse-complement',
                help=(
                    'Include translation of the reverse complement of the '
                    'provided sequence(s) during discovery.'
                ),
            ),
        ] = False,
        enable_gpu: Annotated[ORFaqsCli._enable_gpu_annotation()] = False,
        launch_json: Annotated[ORFaqsCli._launch_json_annotation()] = None,
        ctx: typer.Context = None,
    ):
        ui_kwargs = ORFaqsCli._process_ui_kwargs(**ctx.params)

        if ui_kwargs.get('input_sequence') is None:
            message = (
                '[ERROR] No input sequence or sequence files were specified. '
                'An input sequence or sequence file must be given either '
                'through its positional argument, or in the launch json.'
            )
            _logger.error(message)
            typer.echo(ctx.get_help())
            raise typer.Exit(
                ORFaqsCliExitCodes.ERROR_REQUIRED_INPUT_NOT_SPECIFIED
            )

        ORFaqsProteinDiscoveryCli._discover_proteins(**ui_kwargs)

    @staticmethod
    def _run():
        app()


if __name__ == '__main__':
    multiprocessing.freeze_support()
    ORFaqsProteinDiscoveryCli.run()
