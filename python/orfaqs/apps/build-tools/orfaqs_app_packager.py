#!/usr/bin/env python

import logging
import os
import pathlib
import PyInstaller.__main__

from tqdm import tqdm

from orfaqs.apps.common.cliutils import CliUtil
from orfaqs.apps.common.orfaqsapp import ORFaqsApp

from orfaqs.lib.utils.directoryutils import DirectoryUtils

_logger = logging.getLogger(__name__)

_AVAILABLE_PACKAGES_NAME_PATH_MAP: dict[str, pathlib.Path] = {}


class ORFaqsAppsPackager(ORFaqsApp):
    """ORFaqsAppsPackager"""

    _IGNORE_PATHS = ['build-tools', 'common']

    @staticmethod
    def program_name():
        return 'ORFaqs Apps Packager'

    @staticmethod
    def _ignore_path(file_path: str | os.PathLike):
        file_path = DirectoryUtils.make_path_object(file_path)
        for path in ORFaqsAppsPackager._IGNORE_PATHS:
            if path in str(file_path):
                return True

        return False

    @staticmethod
    def _package_name_path_map() -> dict[str, pathlib.Path]:
        # Retrieve the ORFaqs app directory from this app's path.
        global _AVAILABLE_PACKAGES_NAME_PATH_MAP
        if not _AVAILABLE_PACKAGES_NAME_PATH_MAP:
            orfaqs_app_directory = DirectoryUtils.make_path_object(
                __file__
            ).parent.parent
            file_path_list = DirectoryUtils.glob_files(
                orfaqs_app_directory, pattern='*.py', recursive=True
            )
            for file_path in file_path_list:
                if ORFaqsAppsPackager._ignore_path(file_path):
                    continue
                else:
                    package_name = DirectoryUtils.remove_extension(
                        file_path
                    ).name
                    _AVAILABLE_PACKAGES_NAME_PATH_MAP[package_name] = file_path

        return _AVAILABLE_PACKAGES_NAME_PATH_MAP

    @staticmethod
    def available_packages() -> list[str]:
        return list(ORFaqsAppsPackager._package_name_path_map().keys())

    @staticmethod
    def print_available_packages():
        packages_str = ', '.join(ORFaqsAppsPackager.available_packages())
        print(f'Available Packages:\n\t{packages_str}')

    @staticmethod
    def build_package(package_name: str) -> str:
        app_path = ORFaqsAppsPackager._package_name_path_map().get(
            package_name
        )
        if app_path is None:
            return

        app_path_str = str(app_path)
        build_command_args = [app_path_str, '--paths=./orfaqs']
        PyInstaller.__main__.run(build_command_args)

    @staticmethod
    def build_all_packages():
        # Find and builld all available packages...
        process_list = tqdm(
            ORFaqsAppsPackager.available_packages(),
            desc='Building packages...',
        )
        for package_name in process_list:
            ORFaqsAppsPackager.build_package(package_name)

    @staticmethod
    def cli():
        try:
            arg_descriptor_list = [
                CliUtil.create_new_arg_descriptor(
                    ('-p', '--package_name'),
                    arg_help=('The name of the package to build.'),
                    arg_type=str,
                    default=None,
                ),
                CliUtil.create_new_arg_descriptor(
                    ('--build_all'),
                    arg_help=('Builds all available packages.'),
                    action='store_true',
                ),
                CliUtil.create_new_arg_descriptor(
                    ('--available_packages'),
                    arg_help=('Prints a list of available packages.'),
                    action='store_true',
                ),
            ]

            cli_arg_parser = CliUtil.create_arg_parser(
                arg_descriptor_list,
                program_name=ORFaqsAppsPackager.program_name(),
                description=('Discover proteins from genomic sequences.'),
                epilog='',
            )
            ui_args = CliUtil.parse_args(cli_arg_parser)
            if (
                (ui_args.get('package_name') is not None)
                and (
                    ui_args.get('build_all')
                    or ui_args.get('available_packages')
                )
            ) or (
                ui_args.get('build_all') and ui_args.get('available_packages')
            ):
                message = (
                    '[ERROR] Conflicting arguments. '
                    'Pass only one argument to the CLI'
                )
                print(message)
                cli_arg_parser.print_help()
            elif (
                (ui_args.get('package_name') is None)
                and (not ui_args.get('build_all'))
                and (not ui_args.get('available_packages'))
            ):
                cli_arg_parser.print_help()
                ORFaqsAppsPackager.print_available_packages()
            else:
                ORFaqsAppsPackager.run(**ui_args)
        except SystemExit as error:
            if error.args[0] is True:
                raise

    @staticmethod
    def _run(
        package_name: str = None,
        build_all: bool = False,
        available_packages: bool = False,
    ):
        if package_name is not None:
            ORFaqsAppsPackager.build_package(package_name)
        elif build_all:
            ORFaqsAppsPackager.build_all_packages()
        elif available_packages:
            ORFaqsAppsPackager.print_available_packages()


if __name__ == '__main__':
    ORFaqsAppsPackager.cli()
