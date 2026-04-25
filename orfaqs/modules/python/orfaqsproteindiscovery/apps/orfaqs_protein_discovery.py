#!/usr/bin/env python

import logging
import multiprocessing
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


_logger = logging.getLogger(__name__)
app = ORFaqsCli._new_app_typer()


class ORFaqsProteinDiscoveryCli(ORFaqsCli):
    """ORFaqsProteinDiscovery"""

    @staticmethod
    def program_name():
        return 'ORFaqs Protein Discovery'

    @app.callback(invoke_without_command=True)
    @staticmethod
    def discover_proteins(
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
        if ui_kwargs.get('output_directory') is None:
            ui_kwargs['output_directory'] = (
                ORFaqsProteinDiscoveryCli.default_output_directory()
            )
        ORFaqsProteinDiscoveryApi.discover_proteins(**ui_kwargs)

    @staticmethod
    def _run():
        app()


if __name__ == '__main__':
    multiprocessing.freeze_support()
    ORFaqsProteinDiscoveryCli.run()
