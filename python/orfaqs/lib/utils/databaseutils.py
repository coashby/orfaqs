"""
Database Utils module supporting database connections, queries, and other
common operations.
"""

import csv
import io
import logging
import sqlalchemy
import typing

from abc import ABC, abstractmethod
from pandas.io.sql import SQLTable
from sqlalchemy.exc import (
    ProgrammingError,
)

from sqlalchemy.orm import DeclarativeBase
from sqlalchemy.orm import mapped_column


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


class BaseTable(DeclarativeBase):
    """BaseTable
    Base class for declarative table creation and other declarative mappings.
    """

    pass


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
        query: str,
        connection: sqlalchemy.Connection,
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
    def create_database(
        database: str,
        connection: sqlalchemy.Connection,
    ):
        query = f'CREATE DATABASE {database}'
        connection.execution_options(isolation_level='AUTOCOMMIT')
        SqlAlchemyUtils.run_query(query, connection)
        message = f'[INFO] Database: {database} created.'
        _logger.info(message)
        print(message)

    @staticmethod
    def drop_database(
        database: str,
        connection: sqlalchemy.Connection,
    ):
        query = f'DROP DATABASE IF EXISTS {database} WITH (FORCE)'
        connection.execution_options(isolation_level='AUTOCOMMIT')
        SqlAlchemyUtils.run_query(query, connection)
        message = f'[INFO] Database: {database} dropped.'
        _logger.info(message)
        print(message)

    @staticmethod
    def table_exists(
        table: str,
        engine: sqlalchemy.Engine,
    ) -> bool:
        database_inspector = sqlalchemy.inspect(engine)
        return database_inspector.has_table(table)

    @staticmethod
    def drop_table(
        table: str,
        connection: sqlalchemy.Connection,
    ):
        query = f'DROP TABLE IF EXISTS {table}'
        connection.execution_options(isolation_level='AUTOCOMMIT')
        SqlAlchemyUtils.run_query(query, connection)
        message = f'[INFO] Table: {table} dropped.'
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

    create_column: typing.Callable = mapped_column

    @staticmethod
    def psql_insert_copy(
        table: SQLTable,
        connection: (sqlalchemy.Connection | sqlalchemy.Engine),
        columns: list[str],
        data_iter: typing.Iterable,
    ):
        """
        Inserts data into the specified SQL table.
        ---------
        Arguments
        ---------
        table (SQLTable):
        connection (sqlalchemy.engine.Connection):
        keys (list[str]):
            Column names
        data_iter (typing.Iterable):
        """
        # Retrieve a cursor object from the database API connection.
        database_connection = connection.connection
        with database_connection.cursor() as cursor:
            stream_buffer = io.StringIO()
            csv_writer = csv.writer(stream_buffer)
            csv_writer.writerows(data_iter)
            stream_buffer.seek(0)

            columns_str = ', '.join([f'"{column}"' for column in columns])
            if table.schema:
                table_name = f'{table.schema}.{table.name}'
            else:
                table_name = table.name

            sql = f'COPY {table_name} ({columns_str}) FROM STDIN WITH CSV'
            cursor.copy_expert(sql=sql, file=stream_buffer)


class DatabaseUtils(ABC):
    @staticmethod
    def supported_dialects() -> list[str]:
        return _SUPPORTED_DIALECTS

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
    def database_exists(database) -> bool:
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
    def drop_database(
        database: str,
        **kwargs,
    ):
        pass

    @classmethod
    def table_exists(
        cls,
        table: str,
        connection_options: DatabaseConnectionOptions = None,
    ) -> bool:
        engine = cls.get_engine(connection_options)
        return SqlAlchemyUtils.table_exists(table, engine)

    @staticmethod
    @abstractmethod
    def create_table(
        base_table: DeclarativeBase,
        connection_options: DatabaseConnectionOptions = None,
        **kwargs,
    ):
        pass

    @staticmethod
    @abstractmethod
    def drop_table(
        table: str,
        **kwargs,
    ):
        pass

    @staticmethod
    def run_query(
        query: str,
        connection: sqlalchemy.Connection,
    ) -> sqlalchemy.CursorResult[any]:
        return SqlAlchemyUtils.run_query(query, connection)


class PostgresDatabaseUtils(DatabaseUtils):
    """PostgresDatabaseUtils"""

    @staticmethod
    def default_database() -> str:
        return 'postgres'

    @staticmethod
    def default_username() -> str:
        return 'postgres'

    @staticmethod
    def default_database_connection_options():
        return DatabaseConnectionOptions(
            database=PostgresDatabaseUtils.default_database(),
            username=PostgresDatabaseUtils.default_username(),
            host='localhost',
            port='5432',
        )

    @staticmethod
    def get_engine(
        connection_options: DatabaseConnectionOptions = None,
        **kwargs,
    ) -> sqlalchemy.Engine:
        if connection_options is None:
            connection_options = (
                PostgresDatabaseUtils.default_database_connection_options()
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
        default_database_connection = PostgresDatabaseUtils.connect()
        query = (
            'SELECT EXISTS ('
            f"SELECT 1 FROM pg_database WHERE datname = '{database}')"
        )
        result = SqlAlchemyUtils.run_query(query, default_database_connection)
        default_database_connection.close()
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
            default_database_connection = PostgresDatabaseUtils.connect()
            SqlAlchemyUtils.create_database(
                connection_options.database,
                default_database_connection,
            )
            default_database_connection.close()

    @staticmethod
    def drop_database(database, **kwargs):
        database_connection = PostgresDatabaseUtils.connect(**kwargs)
        SqlAlchemyUtils.drop_database(database, database_connection)
        database_connection.close()

    @staticmethod
    def create_table(
        base_table: DeclarativeBase,
        connection_options: DatabaseConnectionOptions = None,
        **kwargs,
    ):
        database_engine = PostgresDatabaseUtils.get_engine(connection_options)
        base_table.metadata.create_all(database_engine)

    @staticmethod
    def drop_table(
        table: str,
        connection_options: DatabaseConnectionOptions = None,
        **kwargs,
    ):
        database_connection = PostgresDatabaseUtils.connect(connection_options)
        SqlAlchemyUtils.drop_table(table, database_connection)
        database_connection.close()
