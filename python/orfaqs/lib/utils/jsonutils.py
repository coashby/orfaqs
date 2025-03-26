'''
JSON Utils
'''

import enum
import json
import logging
import os

from json.decoder import JSONDecodeError

from orfaqs.lib.utils.directoryutils import DirectoryUtils

_logger = logging.getLogger(__name__)

_JSON_DEFAULT_INDENT = 4
_JSON_FILE_EXTENSION = '.json'


class JsonUtils:
    '''JsonUtils'''

    @staticmethod
    def json_file_extension():
        return _JSON_FILE_EXTENSION

    @staticmethod
    def json_default_indent():
        return _JSON_DEFAULT_INDENT

    @staticmethod
    def make_writable(value):
        if isinstance(value, (str, int, float)):
            return value
        elif isinstance(value, bytes):
            return value.decode()
        elif isinstance(value, enum.Enum):
            return value.name
        elif isinstance(value, dict):
            writable_object = {}
            for key, sub_value in value.items():
                writable_object[key] = JsonUtils.make_writable(sub_value)
            return writable_object
        elif isinstance(value, list):
            writable_object = []
            for sub_value in value:
                writable_object.append(JsonUtils.make_writable(sub_value))
            return writable_object
        else:
            try:
                value = str(value)
            except ValueError as error:
                message = ('[ERROR] Could not find an appropriate conversion'
                           f'for the input value... {error}')
                _logger.error(message)
                return None

        return value

    @staticmethod
    def read_json(json_input: str | os.PathLike,
                  raise_exceptions: bool = True):
        json_contents = None
        if DirectoryUtils.is_file(json_input):
            try:
                with open(json_input, 'r', encoding='utf-8') as i_file:
                    json_contents = json.load(i_file)
            except IOError as error:
                message = (
                    '[ERROR: IOError] '
                    f'An exception occurred while reading from: "{json_input}"'
                )
                _logger.info(message)
                _logger.error(error)
                if raise_exceptions:
                    raise error
            except JSONDecodeError as error:
                message = (
                    '[ERROR: JSONDecodeError] '
                    f'An exception occurred while reading from: "{json_input}"'
                )
                _logger.info(message)
                _logger.error(error)

                if raise_exceptions:
                    raise error

        else:
            try:
                json_contents = json.loads(json_input)

            except JSONDecodeError as error:
                message = (
                    f'[ERROR: JSONDecodeError] {error}'
                )
                _logger.error(message)

                if raise_exceptions:
                    raise error

        return json_contents

    @staticmethod
    def write_json(file_path_str,
                   content,
                   mode='w',
                   indent=None):
        try:
            file_path = DirectoryUtils.make_path_object(file_path_str)
            DirectoryUtils.mkdir_path(file_path.parent)
            with open(file_path_str, mode, encoding='utf-8') as o_file:
                o_file.write(json.dumps(content, indent=indent))
        except IOError as error:
            message = (
                f'An exception occurred while reading from: "{file_path_str}"'
            )
            _logger.info(message)
            _logger.error(error)

    @staticmethod
    def write_json_pretty_print(file_path_str, content, mode='w'):
        JsonUtils.write_json(file_path_str, content, mode,
                             indent=_JSON_DEFAULT_INDENT)

    @staticmethod
    def as_json_string(content, indent=_JSON_DEFAULT_INDENT):
        return json.dumps(content, indent=indent)
