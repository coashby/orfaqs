"""
Database Utils test module.
"""

import logging
import os
import pandas as pd
import pytest
from orfaqs.lib.utils.databaseutils import (
    BaseTable,
    DatabaseConnectionOptions,
    PostgresDatabaseUtils,
    sqlalchemy,
    SqlAlchemyUtils,
)

_logger = logging.getLogger(__name__)


def _database_connection_options(database=None) -> DatabaseConnectionOptions:
    # Check the systems environment variables for database credentials.
    username = os.environ.get('DATABASE_USERNAME')
    password = os.environ.get('DATABASE_PASSWORD')
    host = os.environ.get('DATABASE_HOST')
    port = os.environ.get('DATABASE_PORT')

    return DatabaseConnectionOptions(
        database=database,
        username=username,
        password=password,
        host=host,
        port=port,
    )


def _define_test_table(table: str, id_column_key: str) -> str:
    if table not in BaseTable.metadata.tables:

        class _TestTable(BaseTable):
            """_TestTable"""

            __tablename__ = table
            id = SqlAlchemyUtils.create_column(
                id_column_key,
                sqlalchemy.String,
                primary_key=True,
            )

    return table


class _TestParams:
    def __init__(
        self,
        database_connection_options: DatabaseConnectionOptions = None,
        test_table: str = None,
        test_table_id_key: str = None,
    ):
        self._database_connection_options = database_connection_options
        self._test_table = test_table
        self._test_table_id_key = test_table_id_key

    @property
    def database_connection_options(self) -> DatabaseConnectionOptions:
        return self._database_connection_options

    @property
    def database(self) -> str:
        return self._database_connection_options.database

    @property
    def test_table(self) -> str:
        return self._test_table

    @property
    def test_table_id_key(self) -> str:
        return self._test_table_id_key


@pytest.fixture
def database_resources():
    test_database = 'pytest_database'
    test_table = 'pytest_table'
    test_table_id_key = 'id'
    database_connection_options = _database_connection_options(test_database)
    # Create the database iff it does not exist.
    if not PostgresDatabaseUtils.database_exists(test_database):
        PostgresDatabaseUtils.create_database(database_connection_options)

    test_table = _define_test_table(
        table=test_table,
        id_column_key=test_table_id_key,
    )
    if not PostgresDatabaseUtils.table_exists(
        test_table,
        database_connection_options,
    ):
        PostgresDatabaseUtils.create_table(
            BaseTable,
            database_connection_options,
        )

    yield _TestParams(
        database_connection_options,
        test_table,
    )

    message = '[INFO] Running database_resources test clean up.'
    _logger.info(message)
    if PostgresDatabaseUtils.database_exists(test_database):
        # Remove the database.
        PostgresDatabaseUtils.drop_database(test_database)


class TestPostgresDatabaseUtils:
    """TestPostgresDatabaseUtils"""

    @staticmethod
    def test_create_database(database_resources: _TestParams):
        # Asset that the new database exists.
        assert PostgresDatabaseUtils.database_exists(
            database_resources.database
        )

    @staticmethod
    def test_drop_database(database_resources: _TestParams):
        # Remove the database.
        PostgresDatabaseUtils.drop_database(database_resources.database)

        # Assert that the database has been removed.
        assert not PostgresDatabaseUtils.database_exists(
            database_resources.database
        )

    @staticmethod
    def test_create_table(database_resources: _TestParams):
        # Asset that the new database exists.
        assert PostgresDatabaseUtils.table_exists(
            database_resources.test_table,
            database_resources.database_connection_options,
        )

    @staticmethod
    def test_drop_table(database_resources: _TestParams):
        PostgresDatabaseUtils.drop_table(
            database_resources.test_table,
            database_resources.database_connection_options,
        )

        assert not PostgresDatabaseUtils.table_exists(
            database_resources.test_table,
            database_resources.database_connection_options,
        )

    @staticmethod
    def test_table_insert_dataframe(database_resources: _TestParams):
        # Create the DataFrame object.
        number_ids = int(1e3)
        test_dataframe = pd.DataFrame()
        test_dataframe[database_resources.test_table_id_key] = list(
            range(number_ids)
        )
        database_connection = PostgresDatabaseUtils.connect(
            database_resources.database_connection_options
        )
        test_dataframe.to_sql(
            database_resources.test_table,
            database_connection,
            if_exists='replace',
            index=False,
            method=SqlAlchemyUtils.psql_insert_copy,
        )
        count_query = f'SELECT COUNT(*) FROM {database_resources.test_table}'
        cursor_result = PostgresDatabaseUtils.run_query(
            count_query,
            database_connection,
        )
        assert cursor_result.scalar_one() == number_ids
        database_connection.close()
