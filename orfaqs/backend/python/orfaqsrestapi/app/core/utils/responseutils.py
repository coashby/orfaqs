"""
Response utility module for handling file streaming and downloads.

This module provides utilities for converting files and directories into
streaming responses, including support for large file handling with chunked
transmission and automatic compression of directories.
"""

import mimetypes
import os

from typing import Generator
from fastapi.responses import StreamingResponse

from orfaqs.libs.python.utils.directoryutils import DirectoryUtils


class ResponseUtils:
    LARGE_PATH_SIZE_THRESHOLD_BYTES = 512 * 1024 * 1024
    DEFAULT_STREAM_CHUNK_SIZE_BYTES = 64 * 1024 * 1024

    @staticmethod
    def is_large_path(path: any) -> bool:
        path_size_bytes = DirectoryUtils.path_size(path)
        return path_size_bytes > ResponseUtils.LARGE_PATH_SIZE_THRESHOLD_BYTES

    @staticmethod
    def _file_to_bytes(
        file_path: str | os.PathLike, chunk_size_bytes: int = None
    ) -> Generator[bytes, None, None]:
        if chunk_size_bytes is None:
            chunk_size_bytes = ResponseUtils.DEFAULT_STREAM_CHUNK_SIZE_BYTES
        with open(file_path, 'rb') as i_file:
            while chunk := i_file.read(chunk_size_bytes):
                yield chunk

    @staticmethod
    def stream_response_from_directory_path(
        directory_path: str | os.PathLike,
        download_file_name: str,
    ):
        archived_directory_path = DirectoryUtils.zip_path(directory_path)
        return StreamingResponse(
            content=ResponseUtils._file_to_bytes(archived_directory_path),
            media_type=mimetypes.types_map['.zip'],
            headers={
                'Content-Disposition': f'attachment; filename="{download_file_name}"'
            },
        )
