import os
import pathlib

from datetime import datetime

from orfaqs.lib.python.utils.directoryutils import DirectoryUtils
from orfaqs.lib.python.utils.jsonutils import JsonUtils
from orfaqs.lib.python.utils.platformutils import PlatformUtils


class ORFaqsApiSessionData:
    SESSION_DATA_KEY = 'session_data'
    SESSION_PREVIOUS_SESSION_TIMESTAMP_KEY = 'previous_session_timestamp'
    SESSION_TIMESTAMP_KEY = 'session_timestamp'

    def __init__(
        self,
        session_data: any = None,
        previous_session_timestamp: str | datetime = None,
        session_timestamp: str | datetime = None,
        **kwargs,
    ):
        if session_data is None:
            session_data: dict = {}
        if session_timestamp is None:
            session_timestamp = datetime.now()
        self._previous_session_timestamp = str(previous_session_timestamp)
        self._session_data = session_data
        self._session_timestamp = str(session_timestamp)

    @property
    def session_timestamp(self) -> str:
        return self._session_timestamp

    @property
    def previous_session_timestamp(self) -> str:
        return self._previous_session_timestamp

    @property
    def session_data(self) -> any:
        return self._session_data

    def as_dict(self) -> dict:
        return {
            ORFaqsApiSessionData.SESSION_PREVIOUS_SESSION_TIMESTAMP_KEY: self.previous_session_timestamp,
            ORFaqsApiSessionData.SESSION_TIMESTAMP_KEY: self._session_timestamp,
            ORFaqsApiSessionData.SESSION_DATA_KEY: self._session_data,
        }


class ORFaqsApi:
    @classmethod
    def app_directory(cls) -> pathlib.Path:
        path_str = None
        if PlatformUtils.is_windows():
            path_str = DirectoryUtils.make_path_object(os.getenv('APPDATA'))
        else:
            path_str = pathlib.Path.home()

        return (
            DirectoryUtils.make_path_object(path_str)
            / '.orfaqs'
            / f'{cls.__name__.lower()}'
        )

    @classmethod
    def _create_app_directory(cls):
        DirectoryUtils.mkdir_path(cls.app_directory())

    @classmethod
    def _session_file_name(cls) -> str:
        return f'{cls.__name__.lower()}-session.json'

    @classmethod
    def _session_file_path(cls) -> str:
        return cls.app_directory() / cls._session_file_name()

    @classmethod
    def _load_session_file(cls) -> ORFaqsApiSessionData:
        kwargs: dict = JsonUtils.read_json(
            cls._session_file_path(), raise_exceptions=False
        )
        if kwargs is None:
            return None

        return ORFaqsApiSessionData(**kwargs)

    @classmethod
    def _update_session_file(cls, session_data: dict):
        previous_session_data_object = cls._load_session_file()
        previous_session_timestamp: str = None
        if isinstance(previous_session_data_object, ORFaqsApiSessionData):
            previous_session_timestamp = (
                previous_session_data_object.session_timestamp
            )

        new_session_data = ORFaqsApiSessionData(
            session_data=session_data,
            previous_session_timestamp=previous_session_timestamp,
        )
        cls._create_app_directory()
        JsonUtils.write_json_pretty_print(
            cls._session_file_path(),
            new_session_data.as_dict(),
        )
