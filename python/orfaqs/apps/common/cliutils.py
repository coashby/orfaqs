import textwrap

from argparse import ArgumentParser


class CliUtil:
    _DEFAULT_ARG_HELP_TEXT_WIDTH = 80

    @staticmethod
    def _format_text(text):

        return '\n'.join(textwrap.wrap(
            text,
            width=CliUtil._DEFAULT_ARG_HELP_TEXT_WIDTH,
            break_long_words=False))

    @staticmethod
    def create_new_arg_descriptor(
            name_or_flags,
            arg_help=None,
            action=None,
            arg_type=None,
            default=None):

        arg_descriptor = {}

        if arg_help and isinstance(arg_help, str):
            arg_descriptor['help'] = CliUtil._format_text(arg_help)

        if action and isinstance(action, str):
            arg_descriptor['action'] = action

        if arg_type:
            arg_descriptor['type'] = arg_type

        if default:
            arg_descriptor['default'] = default

        if isinstance(name_or_flags, str):
            name_or_flags = [name_or_flags]
        return (name_or_flags, arg_descriptor)

    @staticmethod
    def create_arg_parser(
            arg_descriptor_list,
            program_name=None,
            description=None,
            epilog=None):
        arg_parser = ArgumentParser(
            prog=program_name,
            description=description,
            epilog=epilog
        )
        for name_or_flags, arg_descriptor in arg_descriptor_list:
            arg_parser.add_argument(
                *name_or_flags, **arg_descriptor)

        return arg_parser

    @staticmethod
    def parse_args(arg_parser) -> dict:
        return vars(arg_parser.parse_args())
