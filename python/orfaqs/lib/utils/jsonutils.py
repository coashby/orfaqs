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
    @staticmethod
    def json_file_extension():
        return _JSON_FILE_EXTENSION

    @staticmethod
    def json_default_indent():
        return _JSON_DEFAULT_INDENT

    @staticmethod
    def make_writable(value: any) -> any:
        if isinstance(value, (str, int, float)):
            return value
        elif isinstance(value, bytes):
            return value.decode()
        elif isinstance(value, enum.Enum):
            return value.name
        elif isinstance(value, dict):
            writable_object = {}
            for key, sub_value in value.items():
                key = JsonUtils.make_writable(key)
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
                message = (
                    '[ERROR] Could not find an appropriate conversion'
                    f'for the input value... {error}'
                )
                _logger.error(message)
                return None

        return value

    @staticmethod
    def read_json_str(json_str: str, raise_exceptions: bool = True) -> any:
        json_contents = None
        try:
            json_contents = json.loads(json_str)

        except JSONDecodeError as error:
            message = f'[ERROR: JSONDecodeError] {error}'
            _logger.error(message)

            if raise_exceptions:
                raise error

        return json_contents

    @staticmethod
    def read_json(
        json_input: str | os.PathLike, raise_exceptions: bool = True
    ) -> any:
        json_contents = None

        def _log_error(error: Exception):
            message = (
                f'[ERROR: {error}] '
                f'An exception occurred while reading from: "{json_input}"'
            )
            _logger.error(message)

        error: Exception = None
        if DirectoryUtils.is_file(json_input):
            try:
                with open(json_input, 'r', encoding='utf-8') as i_file:
                    json_contents = json.load(i_file)
            except IOError as e:
                error = e
            except JSONDecodeError as e:
                error = e
        else:
            try:
                json_contents = JsonUtils.read_json_str(
                    json_input, raise_exceptions
                )
            except TypeError as e:
                error = e

        if error is not None and raise_exceptions:
            _log_error(error)
            raise error

        return json_contents

    @staticmethod
    def write_json(file_path_str, content, mode='w', indent=None):
        try:
            file_path = DirectoryUtils.make_path_object(file_path_str)
            DirectoryUtils.mkdir_path(file_path.parent)
            with open(file_path_str, mode, encoding='utf-8') as o_file:
                o_file.write(JsonUtils.as_json_string(content, indent))
        except IOError as error:
            message = (
                f'An exception occurred while reading from: "{file_path_str}"'
            )
            _logger.info(message)
            _logger.error(error)

    @staticmethod
    def write_json_pretty_print(file_path_str, content, mode='w'):
        JsonUtils.write_json(
            file_path_str, content, mode, indent=_JSON_DEFAULT_INDENT
        )

    @staticmethod
    def as_json_string(content, indent=_JSON_DEFAULT_INDENT):
        json_str = None
        try:
            json_str = json.dumps(content, indent=indent)
        except TypeError:
            json_str = json.dumps(
                JsonUtils.make_writable(content), indent=indent
            )

        return json_str
