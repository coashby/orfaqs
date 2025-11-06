"""
Database Utils module supporting database connections, queries, and other
common operations.
"""

import logging
import sqlalchemy
import typing

from abc import ABC, abstractmethod
from sqlalchemy.exc import (
    ProgrammingError,
)

_logger = logging.getLogger(__name__)


class SupportedDialects:
    POSTGRES = 'postgres'


_DialectOptions = typing.Literal[SupportedDialects.POSTGRES]
_SUPPORTED_DIALECTS = [
    database_type for database_type in _DialectOptions.__args__
]


class DatabaseConnectionOptions:
    """DatabaseConnectionOptions"""

    def __init__(
        self,
        database: str = None,
        username: str = None,
        password: str = None,
        host: str = None,
        port: str = None,
    ):
        self._database = database
        self._username = username
        self._password = password
        self._host = host
        self._port = port

    @property
    def database(self) -> str:
        return self._database

    @property
    def username(self) -> str:
        return self._username

    @property
    def password(self) -> str:
        return self._password

    @property
    def host(self) -> str:
        return self._host

    @property
    def port(self) -> str:
        return self._port

    @database.setter
    def database(self, value: str):
        self._database = value

    @username.setter
    def username(self, value: str):
        self._username = value

    @password.setter
    def password(self, value: str):
        self._password = value

    @host.setter
    def host(self, value: str):
        self._host = value

    @port.setter
    def port(self, value: str):
        self._port = value

    @property
    def kwargs(self) -> dict[str, str]:
        return {
            'database': self._database,
            'username': self._username,
            'password': self._password,
            'host': self._host,
            'port': self._port,
        }


class SqlAlchemyUtils:
    """SqlAlchemyUtils"""

    @staticmethod
    def create_engine(engine_url: str) -> sqlalchemy.Engine:
        try:
            return sqlalchemy.create_engine(engine_url)
        except ProgrammingError as e:
            message = f'[ERROR] Could not connect to PostgresSQL engine.{e}'
            _logger.error(message)
            raise RuntimeError(message) from e

    @staticmethod
    def run_query(
        connection: sqlalchemy.Connection, query: str
    ) -> sqlalchemy.CursorResult[any]:
        try:
            return connection.execute(sqlalchemy.text(query))
        except Exception as e:
            message = (
                '[ERROR] An error occurred while executing the user query.\n'
                f'{e}\n'
                '(debug) ->\n'
                f'\tquery: {query}'
            )
            _logger.error(message)
            raise RuntimeError(message) from e

    @staticmethod
    def create_database(connection: sqlalchemy.Connection, database: str):
        query = f'CREATE DATABASE {database}'
        connection.execution_options(isolation_level='AUTOCOMMIT')
        SqlAlchemyUtils.run_query(connection, query)
        message = f'[INFO] Database: {database} created.'
        _logger.info(message)
        print(message)

    @staticmethod
    def drop_database(connection: sqlalchemy.Connection, database: str):
        query = f'DROP DATABASE IF EXISTS {database}'
        connection.execution_options(isolation_level='AUTOCOMMIT')
        SqlAlchemyUtils.run_query(connection, query)
        message = f'[INFO] Database: {database} dropped.'
        _logger.info(message)
        print(message)

    @staticmethod
    def create_postgres_database_url(
        connection_options: DatabaseConnectionOptions = None,
    ):
        if connection_options is None:
            connection_options = DatabaseConnectionOptions()
        dialect = 'postgresql'
        driver = 'psycopg2'
        return sqlalchemy.URL.create(
            drivername=f'{dialect}+{driver}',
            username=connection_options.username,
            password=connection_options.password,
            host=connection_options.host,
            port=connection_options.port,
            database=connection_options.database,
        )


class DatabaseUtils(ABC):
    @staticmethod
    def supported_dialects() -> list[str]:
        return _SUPPORTED_DIALECTS

    @staticmethod
    @abstractmethod
    def database_exists(database) -> bool:
        pass

    @staticmethod
    @abstractmethod
    def get_engine(
        connection_options: DatabaseConnectionOptions = None,
        **kwargs,
    ):
        pass

    @staticmethod
    @abstractmethod
    def connect(
        connection_options: DatabaseConnectionOptions = None,
        **kwargs,
    ):
        pass

    @staticmethod
    @abstractmethod
    def create_database(
        connection_options: DatabaseConnectionOptions = None,
        **kwargs,
    ):
        pass

    @staticmethod
    @abstractmethod
    def drop_database(database, **kwargs):
        pass

    @staticmethod
    def run_query(
        connection: sqlalchemy.Connection, query: str
    ) -> sqlalchemy.CursorResult[any]:
        return SqlAlchemyUtils.run_query(connection, query)


class PostgresDatabaseUtils(DatabaseUtils):
    """PostgresDatabaseUtils"""

    @staticmethod
    def default_database() -> str:
        return 'postgres'

    @staticmethod
    def default_username() -> str:
        return 'postgres'

    @staticmethod
    def get_engine(
        connection_options: DatabaseConnectionOptions = None,
        **kwargs,
    ) -> sqlalchemy.Engine:
        if connection_options is None:
            connection_options = DatabaseConnectionOptions(
                database=PostgresDatabaseUtils.default_database(),
                username=PostgresDatabaseUtils.default_username(),
                host='localhost',
                port='5432',
            )
        engine_url = SqlAlchemyUtils.create_postgres_database_url(
            connection_options
        )

        return SqlAlchemyUtils.create_engine(engine_url)

    @staticmethod
    def connect(
        connection_options: DatabaseConnectionOptions = None,
        **kwargs,
    ) -> sqlalchemy.Connection:
        # Create the PostgresSQL engine.
        database_engine = PostgresDatabaseUtils.get_engine(connection_options)
        return database_engine.connect()

    @staticmethod
    def database_exists(database) -> bool:
        # Connect to the default database.
        default_connection = PostgresDatabaseUtils.connect()
        query = (
            'SELECT EXISTS ('
            f"SELECT 1 FROM pg_database WHERE datname = '{database}')"
        )
        result = SqlAlchemyUtils.run_query(default_connection, query)

        return result.scalar()

    @staticmethod
    def create_database(
        connection_options: DatabaseConnectionOptions = None,
        **kwargs,
    ):
        if connection_options.username is None:
            connection_options.username = (
                PostgresDatabaseUtils.default_username()
            )

        # Test if this database already exists.
        if PostgresDatabaseUtils.database_exists(connection_options.database):
            message = (
                f'[INFO] Database: {connection_options.database} '
                'already exists.'
            )
            _logger.info(message)
            print(message)
        else:
            # Connect to the default PostgresSQL database. The CREATE DATABASE
            # query must be run within a database connection as itt can not be
            # executed inside transactions (as would be the case from
            # engine.execute()).
            default_connection = PostgresDatabaseUtils.connect()
            SqlAlchemyUtils.create_database(
                default_connection,
                connection_options.database,
            )

    @staticmethod
    def drop_database(database, **kwargs):
        default_connection = PostgresDatabaseUtils.connect()
        SqlAlchemyUtils.drop_database(default_connection, database)
