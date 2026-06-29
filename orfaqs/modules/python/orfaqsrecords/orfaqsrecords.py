"""
ORFaqs protein records module.
"""

import abc
import logging
import pandas as pd

from orfaqs.libs.python.core.enzymes import RNAPolymerase
from orfaqs.libs.python.core.nucleotides import (
    GenomicSequence,
    NucleotideUtils,
    RNASequence,
    StrandType,
)
from orfaqs.libs.python.core.proteins import Protein
from orfaqs.libs.python.core.ribosomes import RNAReadingFrame

from orfaqs.libs.python.utils.jsonutils import JsonUtils

_logger = logging.getLogger(__name__)


class ORFaqsProteinRecordKeys:
    UID_KEY = 'uid'
    PROTEIN_KEY = 'protein'
    PROTEIN_LENGTH_KEY = 'protein_length'


class ORFaqsDiscoveredProteinRecordKeys(ORFaqsProteinRecordKeys):
    """ORFaqsDiscoveredProteinRecordKeys"""

    GENOMIC_SEQUENCE_KEY = 'genomic_sequence'
    ORF_N_TERMINUS_KEY = 'orf_n_terminus'
    PROTEIN_KEY = 'protein'
    PROTEIN_LENGTH_KEY = 'protein_length'
    READING_FRAME_KEY = 'reading_frame'
    SOURCE_UID_KEY = 'source_uid'
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
        orf_sequence_position: int,
        orf_sequence: GenomicSequence | str,
        protein: Protein,
    ):
        self._uid = uid
        self._strand_type = strand_type
        self._reading_frame = reading_frame
        self._orf_sequence_position = orf_sequence_position
        if isinstance(orf_sequence, str):
            orf_sequence = NucleotideUtils.make_sequence_object(
                orf_sequence,
                strand_type=strand_type,
            )
        if isinstance(orf_sequence, RNASequence):
            orf_sequence = RNAPolymerase.reverse_transcribe(orf_sequence)
        self._genomic_sequence = orf_sequence
        self._protein = protein

    @property
    def uid(self) -> str:
        return self._uid

    @property
    def strand_type(self) -> StrandType:
        return self._strand_type

    @property
    def reading_frame(self) -> RNAReadingFrame:
        return self._reading_frame

    @property
    def genomic_sequence_position(self) -> int:
        return self._orf_sequence_position

    @property
    def genomic_sequence(self) -> GenomicSequence:
        return self._genomic_sequence

    @property
    def protein(self) -> Protein:
        return self._protein

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
                ORFaqsDiscoveredProteinRecord.ORF_N_TERMINUS_KEY == record_key
            ):
                record_map[record_key] = self._orf_sequence_position
            elif (
                ORFaqsDiscoveredProteinRecord.GENOMIC_SEQUENCE_KEY
                == record_key
            ):
                record_map[record_key] = self._genomic_sequence.sequence_str
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
            ORFaqsDiscoveredProteinRecord.ORF_N_TERMINUS_KEY,
            ORFaqsDiscoveredProteinRecord.GENOMIC_SEQUENCE_KEY,
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
