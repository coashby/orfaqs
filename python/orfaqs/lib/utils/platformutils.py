import platform


class PlatformUtils:
    """PlatformUtils"""

    @staticmethod
    def is_macos() -> bool:
        return platform.system().lower() == 'darwin'

    @staticmethod
    def is_windows() -> bool:
        return platform.system().lower() == 'windows'

    @staticmethod
    def is_linux() -> bool:
        return platform.system().lower() == 'linux'
