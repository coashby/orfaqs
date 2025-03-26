'''
Directory Utils
'''
import os
import pathlib
import shutil


class DirectoryUtils:
    '''DirectoryUtils'''

    @staticmethod
    def is_directory(path):
        return pathlib.Path(path).is_dir()

    @staticmethod
    def is_file(path):
        return pathlib.Path(path).is_file()

    @staticmethod
    def current_directory_path():
        return pathlib.Path.cwd()

    @staticmethod
    def make_path_object(path):
        if isinstance(path, pathlib.Path):
            return path

        if isinstance(path, list):
            path = '/'.join(path)

        return pathlib.Path(path)

    @staticmethod
    def mkdir_path(path):
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def remove_file_path(path):
        if pathlib.Path(path).is_file():
            os.remove(path)

    @staticmethod
    def path_exists(path):
        return pathlib.Path(path).exists()

    @staticmethod
    def glob_files(path, pattern='*.*', recursive=False):
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
    def append_extension(path, extension):
        path_object = DirectoryUtils.make_path_object(path)
        if extension[0] != '.':
            extension = f'.{extension}'

        path_as_string = str(path_object) + extension
        return pathlib.Path(path_as_string)

    @staticmethod
    def remove_extension(path):
        output_path = DirectoryUtils.make_path_object(path)
        file_name_str = str(output_path.name)
        # Find the last occurrence of '.' for an extension
        for i in range(len(file_name_str)):
            r_index = -(i+1)
            if file_name_str[r_index] == '/':
                break
            if file_name_str[r_index] == '.':
                output_path = output_path.parent / file_name_str[:r_index]
                break

        return output_path
