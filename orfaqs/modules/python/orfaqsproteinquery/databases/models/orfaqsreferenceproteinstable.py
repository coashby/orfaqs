""" """

import logging

from orfaqs.modules.python.orfaqsproteindiscovery.orfaqsproteindiscovery import (
    RNAReadingFrame,
)
from orfaqs.modules.python.orfaqsrecords.orfaqsrecords import (
    ORFaqsProteinRecordKeys,
)

from orfaqs.modules.python.orfaqsproteinquery.databases.orfaqsproteintableutils import (
    ORFaqsProteinTableUtils,
)

from orfaqs.libs.python.utils.databaseutils import (
    BaseTable,
    sqlalchemy,
    SqlAlchemyUtils,
)

_logger = logging.getLogger(__name__)


class ORFaqsReferenceProteinsTableSchema(ORFaqsProteinRecordKeys):
    UID_KEY = 'uid'

    @staticmethod
    def columns() -> list[str]:
        return [
            ORFaqsReferenceProteinsTableSchema.UID_KEY,
            ORFaqsReferenceProteinsTableSchema.PROTEIN_KEY,
            ORFaqsReferenceProteinsTableSchema.PROTEIN_LENGTH_KEY,
        ]


class ORFaqsReferenceProteinsTableFactory:
    """_ORFaqsReferenceProteinsTableFactory"""

    table = 'reference_proteins'

    @staticmethod
    def define_table(table: str = None) -> str:
        """
        Defines the discovered_proteins table when called iff the table
        definition does NOT exists in the base class.
        """
        if table is None:
            table = ORFaqsReferenceProteinsTableFactory.table
        if table in BaseTable.metadata.tables:
            return

        class _ORFaqsReferenceProteinsTable(BaseTable):
            """_ORFaqsReferenceProteinsTable"""

            __tablename__ = table
            UID_MAX_CHAR_LENGTH = 64

            @staticmethod
            def _reading_frame_check_constraint() -> str:
                reading_frame_values = [
                    reading_frame.value for reading_frame in RNAReadingFrame
                ]
                return f'reading_frame IN {tuple(reading_frame_values)}'

            uid = SqlAlchemyUtils.create_column(
                ORFaqsReferenceProteinsTableSchema.UID_KEY,
                sqlalchemy.String(UID_MAX_CHAR_LENGTH),
                primary_key=True,
                comment=ORFaqsProteinTableUtils.uid_comment(),
            )

            protein = SqlAlchemyUtils.create_column(
                ORFaqsReferenceProteinsTableSchema.PROTEIN_KEY,
                sqlalchemy.String,
                nullable=False,
                comment=ORFaqsProteinTableUtils.protein_comment(),
            )

            protein_length = SqlAlchemyUtils.create_column(
                ORFaqsReferenceProteinsTableSchema.PROTEIN_LENGTH_KEY,
                sqlalchemy.Integer,
                nullable=False,
                comment=ORFaqsProteinTableUtils.protein_length_comment(),
            )

        return table
