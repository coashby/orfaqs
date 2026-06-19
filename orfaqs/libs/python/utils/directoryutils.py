"""
Directory Utils
"""

import os
import pathlib
import re
import shutil
import tempfile
import zipfile


class DirectoryUtils:
    """DirectoryUtils"""

    DEFAULT_ZIP_INCREMENTAL_COMPRESSION_CHUNK_SIZE_BYTES = 64 * 1024 * 1024

    @staticmethod
    def is_directory(path) -> bool:
        return pathlib.Path(path).is_dir()

    @staticmethod
    def is_file(path) -> bool:
        return pathlib.Path(path).is_file()

    @staticmethod
    def has_suffix(path) -> bool:
        return len(pathlib.Path(path).suffix) > 0

    @staticmethod
    def current_directory_path() -> pathlib.Path:
        return pathlib.Path.cwd()

    @staticmethod
    def make_path_object(path) -> pathlib.Path:
        if isinstance(path, pathlib.Path):
            return path

        if isinstance(path, list):
            path = '/'.join(path)
        elif isinstance(path, tempfile.TemporaryDirectory):
            path = path.name

        return pathlib.Path(path)

    @staticmethod
    def mkdir_path(path):
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def remove_file_path(path):
        if pathlib.Path(path).is_file():
            os.remove(path)

    @staticmethod
    def path_exists(path) -> bool:
        return pathlib.Path(path).exists()

    @staticmethod
    def glob_files(
        path,
        pattern: str = '*.*',
        recursive: bool = False,
    ) -> list[pathlib.Path]:
        file_list = []
        if recursive:
            file_list = list(
                DirectoryUtils.make_path_object(path).rglob(pattern)
            )
        else:
            file_list = list(
                DirectoryUtils.make_path_object(path).glob(pattern)
            )

        return file_list

    @staticmethod
    def move(source_path, destination_path):
        shutil.move(source_path, destination_path)

    @staticmethod
    def copy(source_path, destination_path):
        shutil.copy(source_path, destination_path)

    @staticmethod
    def append_extension(path, extension) -> pathlib.Path:
        path_object = DirectoryUtils.make_path_object(path)
        if extension[0] != '.':
            extension = f'.{extension}'

        path_as_string = str(path_object) + extension
        return pathlib.Path(path_as_string)

    @staticmethod
    def remove_extension(path) -> pathlib.Path:
        output_path = DirectoryUtils.make_path_object(path)
        file_name_str = str(output_path.name)
        # Find the last occurrence of '.' for an extension
        for i in range(len(file_name_str)):
            r_index = -(i + 1)
            if file_name_str[r_index] == '/':
                break
            if file_name_str[r_index] == '.':
                output_path = output_path.parent / file_name_str[:r_index]
                break

        return output_path

    @staticmethod
    def sanitize_path(path) -> pathlib.Path:
        path = pathlib.Path(path)

        def _sanitize_name(name: str) -> str:
            return re.sub(r'[<>:"/\\|?*]', '', name)

        clean_path = pathlib.Path(path.root)
        path_parts = path.parts
        if len(path.root) > 0:
            path_parts = path_parts[1:]
        for sub_path_name in path_parts:
            clean_path = clean_path / _sanitize_name(sub_path_name)

        return clean_path

    @staticmethod
    def save_as_file(
        file_path: str | pathlib.Path,
        file_contents: any,
        binary_data: bool = False,
    ) -> pathlib.Path:
        file_path = DirectoryUtils.make_path_object(file_path)
        write_mode = 'wb' if binary_data else 'w'
        with open(file_path, write_mode) as o_file:
            shutil.copyfileobj(file_contents, o_file)

        return file_path

    @staticmethod
    def temporary_directory() -> tempfile.TemporaryDirectory:
        return tempfile.TemporaryDirectory()

    @staticmethod
    def file_stats(file_path: any) -> os.stat_result:
        file_path: pathlib.Path = DirectoryUtils.make_path_object(file_path)

        return file_path.stat()

    @staticmethod
    def file_size(file: any) -> int:
        return DirectoryUtils.file_stats(file).st_size

    @staticmethod
    def path_size(path: str | pathlib.Path) -> int:
        path = DirectoryUtils.make_path_object(path)
        if not path.exists():
            return None

        file_list: list[pathlib.Path] = []
        if DirectoryUtils.is_directory(path):
            file_list = DirectoryUtils.glob_files(path, recursive=True)
        else:
            file_list.append(path)

        total_path_size: int = 0
        for file_path in file_list:
            total_path_size += DirectoryUtils.file_size(file_path)

        return total_path_size

    @staticmethod
    def cleanup_temporary_directory(
        temporary_directory: tempfile.TemporaryDirectory,
    ):
        temporary_directory.cleanup()

    @staticmethod
    def zip_path(
        source_path: str | pathlib.Path,
        destination_path: str | pathlib.Path = None,
        chunk_size_bytes: int = None,
    ) -> pathlib.Path:
        if not DirectoryUtils.path_exists(source_path):
            return None
        source_path = DirectoryUtils.make_path_object(source_path)
        if destination_path is None:
            destination_path = source_path.with_suffix('.zip')

        source_path = DirectoryUtils.make_path_object(source_path)
        file_list: list[pathlib.Path] = []
        folder_path = None
        if DirectoryUtils.is_directory(source_path):
            file_list = DirectoryUtils.glob_files(source_path, recursive=True)
            folder_path = source_path
        else:
            file_list.append(source_path)
            folder_path = source_path.parent

        if chunk_size_bytes is None:
            chunk_size_bytes = DirectoryUtils.DEFAULT_ZIP_INCREMENTAL_COMPRESSION_CHUNK_SIZE_BYTES

        with zipfile.ZipFile(
            destination_path, 'w', zipfile.ZIP_DEFLATED
        ) as zipfile_archive:
            for file_path in file_list:
                # Determine the relative path of the file.
                relative_path = os.path.relpath(file_path, folder_path)
                with zipfile_archive.open(
                    relative_path, mode='w'
                ) as o_archive_stream:
                    with open(file_path, 'rb') as i_file_stream:
                        # Use the specified chunk size to stream the archive.
                        while chunk := i_file_stream.read(chunk_size_bytes):
                            o_archive_stream.write(chunk)

        return destination_path
