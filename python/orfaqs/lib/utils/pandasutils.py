"""
Pandas utility module.
"""

import logging
import os
import pandas as pd
import pathlib
import typing

from orfaqs.lib.utils.directoryutils import DirectoryUtils

_logger = logging.getLogger(__name__)


class DataFrameExportFormat:
    CSV = 'csv'
    JSON = 'json'
    XLSX = 'xlsx'


DataFrameExportFormatOptions = typing.Literal[
    DataFrameExportFormat.CSV,
    DataFrameExportFormat.JSON,
    DataFrameExportFormat.XLSX,
]
_AVAILABLE_EXPORT_FORMATS: list[str] = [
    format for format in DataFrameExportFormatOptions.__args__
]


class PandasUtils:
    @staticmethod
    def available_dataframe_export_formats() -> list[str]:
        return _AVAILABLE_EXPORT_FORMATS

    @staticmethod
    def read_file_as_dataframe(file_path: (str | os.PathLike)) -> pd.DataFrame:
        results_dataframe: pd.DataFrame = None
        file_type = DirectoryUtils.make_path_object(file_path).suffix.lower()
        if DataFrameExportFormat.CSV in file_type:
            results_dataframe = pd.read_csv(file_path, index_col=False)
        elif DataFrameExportFormat.JSON in file_type:
            results_dataframe = pd.read_json(file_path)
        elif DataFrameExportFormat.XLSX in file_type:
            results_dataframe = pd.read_excel(file_path, index_col=False)
        else:
            message = (
                '[ERROR] The input file type is not currently supported.\n'
                '(debug) ->\n'
                f'\tfile_path: {file_path}\n'
                f'\tfile_type: {file_type}'
            )
            _logger.error(message)
            raise ValueError(message)

        results_dataframe.reset_index(inplace=True)
        return results_dataframe

    @staticmethod
    def export_dataframe(
        file_path: (str | os.PathLike),
        dataframe: pd.DataFrame,
        export_format: DataFrameExportFormatOptions,
        include_index: bool = False,
        index_label: str = None,
    ) -> pathlib.Path:
        file_path = DirectoryUtils.make_path_object(file_path)
        if include_index and index_label is None:
            index_label = 'index'
        if DataFrameExportFormat.CSV in export_format:
            file_path = file_path.with_suffix(f'.{DataFrameExportFormat.CSV}')
            dataframe.to_csv(
                file_path,
                index_label=index_label,
                index=include_index,
            )
        elif DataFrameExportFormat.JSON in export_format:
            file_path = file_path.with_suffix(f'.{DataFrameExportFormat.JSON}')
            dataframe.to_json(file_path)
        elif DataFrameExportFormat.XLSX in export_format:
            file_path = file_path.with_suffix(f'.{DataFrameExportFormat.XLSX}')
            dataframe.to_excel(
                file_path,
                index_label=index_label,
                index=include_index,
            )
        return file_path
