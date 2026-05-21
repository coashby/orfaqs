""" """

import logging

from orfaqs.modules.python.orfaqsproteinquery.databases.orfaqsproteintableutils import (
    ORFaqsProteinTableUtils,
)
from orfaqs.modules.python.orfaqsproteindiscovery.orfaqsproteindiscovery import (
    RNAReadingFrame,
)
from orfaqs.modules.python.orfaqsrecords.orfaqsrecords import (
    ORFaqsDiscoveredProteinRecordKeys,
)

from orfaqs.libs.python.utils.databaseutils import (
    BaseTable,
    sqlalchemy,
    SqlAlchemyUtils,
)

_logger = logging.getLogger(__name__)


class ORFaqsDiscoveredProteinsTableSchema(ORFaqsDiscoveredProteinRecordKeys):
    UID_KEY = 'uid'

    @staticmethod
    def columns() -> list[str]:
        return [
            ORFaqsDiscoveredProteinsTableSchema.UID_KEY,
            ORFaqsDiscoveredProteinsTableSchema.SOURCE_UID_KEY,
            ORFaqsDiscoveredProteinsTableSchema.STRAND_TYPE_KEY,
            ORFaqsDiscoveredProteinsTableSchema.READING_FRAME_KEY,
            ORFaqsDiscoveredProteinsTableSchema.GENOMIC_SEQUENCE_POSITION_KEY,
            ORFaqsDiscoveredProteinsTableSchema.PROTEIN_KEY,
            ORFaqsDiscoveredProteinsTableSchema.GENOMIC_SEQUENCE_KEY,
            ORFaqsDiscoveredProteinsTableSchema.PROTEIN_LENGTH_KEY,
        ]


class ORFaqsDiscoveredProteinsTableFactory:
    """ORFaqsDiscoveredProteinsTableUtil"""

    table = 'discovered_proteins'

    @staticmethod
    def define_table(table: str = None) -> str:
        """
        Defines the discovered_proteins table when called iff the table
        definition does NOT exists in the base class.
        """
        if table is None:
            table = ORFaqsDiscoveredProteinsTableFactory.table
        if table in BaseTable.metadata.tables:
            return

        class ORFaqsDiscoveredProteinsTable(BaseTable):
            """ORFaqsDiscoveredProteinTable"""

            __tablename__ = table
            UID_MAX_CHAR_LENGTH = 64
            SOURCE_UID_MAX_CHAR_LENGTH = 64
            STRAND_TYPE_MAX_CHAR_LENGTH = 32

            @staticmethod
            def _reading_frame_check_constraint() -> str:
                reading_frame_values = [
                    reading_frame.value for reading_frame in RNAReadingFrame
                ]
                return f'reading_frame IN {tuple(reading_frame_values)}'

            uid = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinsTableSchema.UID_KEY,
                sqlalchemy.String(UID_MAX_CHAR_LENGTH),
                primary_key=True,
                comment=ORFaqsProteinTableUtils.uid_comment(),
            )

            source_uid = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinsTableSchema.SOURCE_UID_KEY,
                sqlalchemy.String(SOURCE_UID_MAX_CHAR_LENGTH),
                comment=ORFaqsProteinTableUtils.source_uid_comment(),
            )

            strand_type = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinsTableSchema.STRAND_TYPE_KEY,
                sqlalchemy.String(STRAND_TYPE_MAX_CHAR_LENGTH),
                nullable=False,
                comment=ORFaqsProteinTableUtils.strand_type_comment(),
            )

            reading_frame = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinsTableSchema.READING_FRAME_KEY,
                sqlalchemy.Integer,
                nullable=False,
                comment=ORFaqsProteinTableUtils.reading_frame_comment(),
            )

            genomic_sequence_position = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinsTableSchema.GENOMIC_SEQUENCE_POSITION_KEY,
                sqlalchemy.Integer,
                sqlalchemy.CheckConstraint(_reading_frame_check_constraint()),
                nullable=False,
                comment=ORFaqsProteinTableUtils.genomic_sequence_position_comment(),
            )

            protein = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinsTableSchema.PROTEIN_KEY,
                sqlalchemy.String,
                nullable=False,
                comment=ORFaqsProteinTableUtils.protein_comment(),
            )

            genomic_sequence = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinsTableSchema.GENOMIC_SEQUENCE_KEY,
                sqlalchemy.String,
                nullable=False,
                comment=ORFaqsProteinTableUtils.genomic_sequence_comment(),
            )

            protein_length = SqlAlchemyUtils.create_column(
                ORFaqsDiscoveredProteinsTableSchema.PROTEIN_LENGTH_KEY,
                sqlalchemy.Integer,
                nullable=False,
                comment=ORFaqsProteinTableUtils.protein_length_comment(),
            )

        return table
