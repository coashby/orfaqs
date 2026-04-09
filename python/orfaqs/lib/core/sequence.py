"""
Sequence module provides abstractions for molecular sequences.
"""

import abc
import logging

_logger = logging.getLogger(__name__)


class Sequence(abc.ABC):
    def __init__(
        self,
        sequence: (str | list),
        name: str,
        log_errors: bool = True,
        raise_errors: bool = True,
    ):
        self._sequence_str: str = ''
        self._name = name
        if isinstance(sequence, list):
            sequence = ''.join(list(map(str, sequence)))

        if isinstance(sequence, str):
            sequence = sequence.upper()

        available_symbols = [
            symbol.upper() for symbol in self.available_symbols()
        ]
        for symbol in sequence:
            symbol = str(symbol)
            if symbol not in available_symbols:
                message = None
                if log_errors:
                    message = (
                        f'[ERROR] "{symbol}" is not a valid symbol for '
                        f'{self.__class__.__name__} objects.'
                    )
                    _logger.error(message)
                if raise_errors:
                    raise ValueError(message)

        self._sequence_str = sequence

    def __contains__(self, region: any) -> bool:
        if isinstance(region, str):
            return region in self._sequence_str
        elif isinstance(region, self.__class__):
            return region._sequence_str in self._sequence_str

        return False

    def __str__(self) -> str:
        return self._sequence_str

    def __len__(self) -> int:
        return len(self._sequence_str)

    def __eq__(self, rhs: any):
        rhs_sequence_str: str = None
        if isinstance(rhs, Sequence):
            rhs_sequence_str = rhs._sequence_str
        elif isinstance(rhs, str):
            rhs_sequence_str = rhs
        return self._sequence_str == rhs_sequence_str

    @property
    def sequence_str(self) -> str:
        return self._sequence_str

    @property
    def sequence_length(self) -> int:
        return len(self._sequence_str)

    @classmethod
    @abc.abstractmethod
    def available_symbols(cls) -> list[str]:
        pass


class SequenceUtils:
    @staticmethod
    def replace_symbols(
        sequence: (str | Sequence), symbol: str, replace_with_symbol: str
    ) -> str:
        sequence_str = sequence
        if isinstance(sequence, Sequence):
            sequence_str = sequence.sequence_str

        return sequence_str.upper().replace(
            symbol.upper(), replace_with_symbol.upper()
        )
