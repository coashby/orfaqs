"""
Database Utils test module.
"""

import os
import pytest
from orfaqs.lib.utils.databaseutils import (
    DatabaseConnectionOptions,
    PostgresDatabaseUtils,
)


@pytest.fixture
def database_connection_options():
    # Check the systems environment variables for database credentials.
    database = os.environ.get('DATABASE_NAME')
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


class TestPostgresDatabaseUtils:
    @staticmethod
    def test_create_database(
        database_connection_options: DatabaseConnectionOptions,
    ):
        test_database = 'test_database'
        database_connection_options.database = test_database
        PostgresDatabaseUtils.create_database(database_connection_options)

        # Asset that the new database exists.
        assert PostgresDatabaseUtils.database_exists(test_database)

        # Remove the database.
        PostgresDatabaseUtils.drop_database(test_database)

        # Assert that the database has been removed.
        assert not PostgresDatabaseUtils.database_exists(test_database)
