"""
ORFaqs protein records module.
"""

import abc
import logging
import pandas as pd

from orfaqs.lib.core.nucleotides import StrandType
from orfaqs.lib.core.proteins import Protein
from orfaqs.lib.core.ribosomes import RNAReadingFrame

from orfaqs.lib.utils.jsonutils import JsonUtils

_logger = logging.getLogger(__name__)


class ORFaqsProteinRecordKeys:
    UID_KEY = 'uid'
    PROTEIN_KEY = 'protein'
    PROTEIN_LENGTH_KEY = 'protein_length'


class ORFaqsDiscoveredProteinRecordKeys(ORFaqsProteinRecordKeys):
    """ORFaqsDiscoveredProteinRecordKeys"""

    SOURCE_UID_KEY = 'source_uid'
    PROTEIN_KEY = 'protein'
    PROTEIN_LENGTH_KEY = 'protein_length'
    READING_FRAME_KEY = 'reading_frame'
    DNA_SEQUENCE_POSITION_KEY = 'dna_sequence_position'
    STRAND_TYPE_KEY = 'strand_type'


class ORFaqsRecord(abc.ABC):
    @property
    @abc.abstractmethod
    def record_dict(self) -> dict[str, any]:
        pass

    @property
    def condensed_record_json_str(self) -> str:
        record_json_str = JsonUtils.as_json_string(
            JsonUtils.make_writable(self.record_dict), indent=None
        )
        record_json_str = record_json_str.strip()
        record_json_str = record_json_str.replace('\n', '')

        return record_json_str

    def pretty_print_record(self):
        record_json_str = JsonUtils.as_json_string(
            JsonUtils.make_writable(self.record_dict)
        )

        print(record_json_str)

    @staticmethod
    @abc.abstractmethod
    def keys() -> list[str]:
        pass


class ORFaqsRecordUtils:
    """ORFaqsRecordUtils"""

    @staticmethod
    def orfaqs_records_to_dataframe(
        records: ORFaqsRecord | list[ORFaqsRecord],
    ) -> pd.DataFrame:
        if not isinstance(records, list):
            records = [records]

        record_dicts: list[dict] = []
        for record in records:
            record_dicts.append(record.record_dict)

        return pd.DataFrame(record_dicts)


class ORFaqsDiscoveredProteinRecord(
    ORFaqsDiscoveredProteinRecordKeys, ORFaqsRecord
):
    """ORFaqsDiscoveredProteinRecord"""

    def __init__(
        self,
        uid: str,
        strand_type: StrandType,
        reading_frame: RNAReadingFrame,
        rna_sequence_position: int,
        protein: Protein,
    ):
        self._uid = uid
        self._strand_type = strand_type
        self._reading_frame = reading_frame
        self._rna_sequence_position = rna_sequence_position
        self._protein = protein

    @property
    def record_dict(self) -> dict[str, any]:
        record_map: dict[str, any] = {}
        for record_key in ORFaqsDiscoveredProteinRecord.keys():
            if ORFaqsDiscoveredProteinRecord.SOURCE_UID_KEY == record_key:
                record_map[record_key] = self._uid
            elif ORFaqsDiscoveredProteinRecord.STRAND_TYPE_KEY == record_key:
                record_map[record_key] = str(self._strand_type)
            elif ORFaqsDiscoveredProteinRecord.READING_FRAME_KEY == record_key:
                record_map[record_key] = self._reading_frame.value
            elif (
                ORFaqsDiscoveredProteinRecord.DNA_SEQUENCE_POSITION_KEY
                == record_key
            ):
                record_map[record_key] = self._rna_sequence_position
            elif ORFaqsDiscoveredProteinRecord.PROTEIN_KEY == record_key:
                record_map[record_key] = self._protein
            elif (
                ORFaqsDiscoveredProteinRecord.PROTEIN_LENGTH_KEY == record_key
            ):
                record_map[record_key] = self._protein.sequence_length
            else:
                message = (
                    '[ERROR] Missing record key assignment.\n'
                    '(debug) ->\n'
                    f'\trecord_key: {record_key} (not assigned)'
                )
                _logger.error(message)
                raise RuntimeError(message)

        return record_map

    @staticmethod
    def keys() -> list[str]:
        return [
            ORFaqsDiscoveredProteinRecord.SOURCE_UID_KEY,
            ORFaqsDiscoveredProteinRecord.STRAND_TYPE_KEY,
            ORFaqsDiscoveredProteinRecord.READING_FRAME_KEY,
            ORFaqsDiscoveredProteinRecord.DNA_SEQUENCE_POSITION_KEY,
            ORFaqsDiscoveredProteinRecord.PROTEIN_KEY,
            ORFaqsDiscoveredProteinRecord.PROTEIN_LENGTH_KEY,
        ]


class ORFaqsProteinRecord(ORFaqsProteinRecordKeys):
    """ORFaqsProteinRecord"""

    def __init__(
        self,
        uid: str,
        protein: Protein,
    ):
        self._uid = uid
        self._protein = protein
        self._protein_length = protein.sequence_length

    @property
    def record_dict(self) -> dict[str, any]:
        record_map: dict[str, any] = {}
        for record_key in ORFaqsProteinRecord.keys():
            if ORFaqsProteinRecord.UID_KEY == record_key:
                record_map[record_key] = self._uid
            elif ORFaqsProteinRecord.PROTEIN_KEY == record_key:
                record_map[record_key] = self._protein
            elif ORFaqsProteinRecord.PROTEIN_LENGTH_KEY == record_key:
                record_map[record_key] = self._protein.sequence_length
            else:
                message = (
                    '[ERROR] Missing record key assignment.\n'
                    '(debug) ->\n'
                    f'\trecord_key: {record_key} (not assigned)'
                )
                _logger.error(message)
                raise RuntimeError(message)

        return record_map

    @staticmethod
    def keys() -> list[str]:
        return [
            ORFaqsProteinRecord.UID_KEY,
            ORFaqsProteinRecord.PROTEIN_KEY,
            ORFaqsProteinRecord.PROTEIN_LENGTH_KEY,
        ]
