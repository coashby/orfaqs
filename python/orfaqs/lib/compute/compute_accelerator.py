import abc
import logging
import os

from numpy.typing import ArrayLike

_logger = logging.getLogger(__name__)


class ComputeAccelerator(abc.ABC):
    """ComputeAccelerator"""

    def __init__(
        self,
        device_name: str = None,
        kernel_source_file_path: (str | os.PathLike) = None,
        **kwargs,
    ):
        """Initializes available compute devices and/or kernels.

        Each ComputeAccelerator manages a single compute device. If no compute
        device is specified, the default compute device is used.
        ---------
        Arguments
        ---------
        device_name (str):
            The name of the compute device to use for computation.

        kernel_source_file_path (str | os.PathLike):
            The file path of the kernel/shader to compile.
        """

        self._compute_device_map: dict[str, any] = {}
        self._compute_device = None
        self._kernel_binary = None
        self._kernel_function_map: dict[str, any] = {}
        self._command_queue: any = None
        self._command_buffer: any = None
        self._number_available_threads: int = 0
        self.init(
            device_name=device_name,
            kernel_source_file_path=kernel_source_file_path,
            **kwargs,
        )

    @property
    def number_available_threads(self) -> int:
        return self._number_available_threads

    @property
    def available_functions(self) -> list[str]:
        return list(self._kernel_function_map.keys())

    @abc.abstractmethod
    def _load_kernel_str(self, kernel_source_str: str) -> any:
        pass

    @abc.abstractmethod
    def _load_kernel_function_map(self, kernel_binary):
        pass

    def _set_device(self, device_name: str):
        if self._kernel_binary is not None:
            message = (
                '[ERROR] A kernel binary already exists. Release the current '
                'kernel before attempting change compute devices.'
            )
            _logger.error(message)
            raise RuntimeError(message)
        if device_name not in self._compute_device_map:
            message = (
                f'[ERROR] The device: "{device_name}" does not exist on '
                'the current platform.'
            )
            raise ValueError(message)

        self._compute_device = self._compute_device_map[device_name]

    @abc.abstractmethod
    def _release(self):
        pass

    @abc.abstractmethod
    def init(
        self,
        device_name: str = None,
        kernel_source_file_path: (str | os.PathLike) = None,
        **kwargs,
    ):
        pass

    @staticmethod
    @abc.abstractmethod
    def platform_supported() -> bool:
        pass

    @staticmethod
    @abc.abstractmethod
    def available_devices() -> list[str]:
        pass

    @staticmethod
    @abc.abstractmethod
    def max_number_threads_on_device(device_name: str = None):
        pass

    def release(self):
        self._release()
        self._kernel_binary = None
        self._compute_device_map = {}
        self._kernel_function_map = {}

    def load_kernel(self, file_path: (str | os.PathLike)):
        kernel_source_str = open(file_path, encoding='utf-8').read()
        self.load_kernel_str(kernel_source_str)

    def load_kernel_str(self, kernel_source_str: str):
        if self._kernel_binary is not None:
            message = (
                '[ERROR] A kernel binary already exists. Release the current '
                'kernel before attempting to load a new kernel.'
            )
            _logger.error(message)
            raise RuntimeError(message)

        self._kernel_binary = self._load_kernel_str(kernel_source_str)
        if self._kernel_binary is None:
            message = (
                '[ERROR] Failed to compile the given kernel. This may '
                'be the result of syntax errors in the kernel source code.'
            )
            _logger.error(message)
            raise ValueError(message)

        self._load_kernel_function_map(self._kernel_binary)

    @abc.abstractmethod
    def set_arg(
        self,
        function_name: str,
        arg_index: int,
        arg: any,
        dtype: type = None,
    ):
        pass

    @abc.abstractmethod
    def execute(self, blocking: bool = True, **kwargs):
        pass

    @abc.abstractmethod
    def get_buffer(self, **kwargs) -> ArrayLike:
        pass
