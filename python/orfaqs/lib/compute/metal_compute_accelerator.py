import logging
import Metal
import numpy as np
import os

from numpy.typing import ArrayLike

from orfaqs.lib.compute.compute_accelerator import ComputeAccelerator

_logger = logging.getLogger(__name__)


class MetalComputeAccelerator(ComputeAccelerator):
    """MetalComputeAccelerator"""

    def __init__(
        self,
        device_name: str = None,
        kernel_source_file_path: (str | os.PathLike) = None,
    ):
        self._command_queue = None
        self._command_buffer = None
        self._current_pipeline_state = None
        self._compute_command_encoder = None
        self._kernel_arg_buffers = None
        super().__init__(
            device_name=device_name,
            kernel_source_file_path=kernel_source_file_path,
        )

    @staticmethod
    def _get_max_number_threads_from_device(compute_device) -> int:
        device_dimensions = compute_device.maxThreadsPerThreadgroup()
        return max(
            [
                device_dimensions.width,
                device_dimensions.height,
                device_dimensions.depth,
            ]
        )

    def _load_kernel_str(self, kernel_source_str: str) -> any:
        return self._compute_device.newLibraryWithSource_options_error_(
            kernel_source_str, None, None
        )[0]

    def _load_kernel_function_map(self, kernel_binary):
        function_names = kernel_binary.functionNames()
        for function_name in function_names:
            self._kernel_function_map[function_name] = (
                kernel_binary.newFunctionWithName_(function_name)
            )

    def _release(self):
        self._compute_device = None

    def init(self, device_name=None, kernel_source_file_path=None, **kwargs):
        available_devices = self.available_devices()
        for device in available_devices:
            self._compute_device_map[device.name()] = device

        if device_name is not None:
            self._set_device(device_name)
        else:
            self._compute_device = Metal.MTLCreateSystemDefaultDevice()

        if self._compute_device is not None:
            # Initialize kernel execution resources.
            self._command_queue = self._compute_device.newCommandQueue()
            self._command_buffer = self._command_queue.commandBuffer()
            self._compute_command_encoder = (
                self._command_buffer.computeCommandEncoder()
            )

        if kernel_source_file_path is not None:
            self.load_kernel(kernel_source_file_path)

        self._number_available_threads = (
            MetalComputeAccelerator._get_max_number_threads_from_device(
                self._compute_device
            )
        )

    @staticmethod
    def platform_supported() -> bool:
        return Metal.MTLCreateSystemDefaultDevice() is not None

    @staticmethod
    def available_devices() -> list[str]:
        return Metal.MTLCopyAllDevices()

    @staticmethod
    def max_number_threads_on_device(device_name: str = None) -> int:
        compute_device = None
        if device_name is None:
            compute_device = Metal.MTLCreateSystemDefaultDevice()
        else:
            for device in MetalComputeAccelerator.available_devices():
                if device.name == device_name:
                    compute_device = device
                    break

        if compute_device is None:
            message = '[ERROR] A valid device was not initialized.'
            if device_name is not None:
                message += f' The device name: "{device_name}" does not exist.'
            _logger.error(message)
            raise ValueError(message)

        return MetalComputeAccelerator._get_max_number_threads_from_device(
            compute_device
        )

    def _configure_compute_resources(self, function_name: str):
        kernel_function = self._kernel_function_map.get(function_name)
        if kernel_function is None:
            message = (
                f'[ERROR] The requested function name: "{function_name}" is '
                'not in any compiled shader libraries.'
            )
            _logger.error(message)
            return ValueError(message)
        # Create the device buffers for the selected kernel function.
        # Create a new pipeline state object to run the requested kernel.
        self._current_pipeline_state = (
            self._compute_device.newComputePipelineStateWithFunction_error_(
                kernel_function,
                None,
            )[0]
        )
        # Set the command encoder's pipeline state to run the requested kernel
        # function.
        self._compute_command_encoder.setComputePipelineState_(
            self._current_pipeline_state
        )

    def set_arg(
        self,
        function_name: str,
        arg_index: int,
        arg: any,
        dtype: type = np.float32,
    ):
        # Configure the compute resources for the given kernel function
        self._configure_compute_resources(function_name)

        # Set all buffer args in order
        buffer_offset = 0
        if self._kernel_arg_buffers is None:
            self._kernel_arg_buffers: list[any] = []
        arg = np.asarray(arg, dtype=dtype)
        buffer_arg = self._compute_device.newBufferWithBytes_length_options_(
            arg,
            arg.nbytes,
            Metal.MTLResourceStorageModeShared,
        )
        self._compute_command_encoder.setBuffer_offset_atIndex_(
            buffer_arg,
            buffer_offset,
            arg_index,
        )
        self._kernel_arg_buffers.insert(arg_index, buffer_arg)

    def set_args(self, function_name: str, *args):
        # Configure the compute resources for the given kernel function
        self._configure_compute_resources(function_name)

        # Set all buffer args in order
        buffer_offset = 0
        self._kernel_arg_buffers: list[any] = []
        for arg_index, arg in enumerate(args):
            # Ensure the input arg is an np.ndarray.
            arg = np.asarray(arg, dtype=np.float32)
            buffer_arg = (
                self._compute_device.newBufferWithBytes_length_options_(
                    arg,
                    arg.nbytes,
                    Metal.MTLResourceStorageModeShared,
                )
            )
            self._compute_command_encoder.setBuffer_offset_atIndex_(
                buffer_arg,
                buffer_offset,
                arg_index,
            )
            self._kernel_arg_buffers.append(buffer_arg)

    def execute(
        self,
        blocking: bool = True,
        threads_per_group: tuple[int] = None,
        number_of_groups: tuple[int] = None,
        **kwargs,
    ):
        if threads_per_group is None:
            threads_per_group = (
                self._current_pipeline_state.maxTotalThreadsPerThreadgroup(),
                1,
                1,
            )
        if number_of_groups is None:
            number_of_groups = (
                self._current_pipeline_state.threadExecutionWidth(),
                1,
                1,
            )
        threads_per_group = Metal.MTLSizeMake(*threads_per_group)
        number_of_groups = Metal.MTLSizeMake(*number_of_groups)

        self._compute_command_encoder.dispatchThreads_threadsPerThreadgroup_(
            threads_per_group,
            number_of_groups,
        )
        self._compute_command_encoder.endEncoding()
        self._command_buffer.commit()

        if blocking:
            self._command_buffer.waitUntilCompleted()

    def get_buffer(
        self,
        buffer_index: int,
        dtype=np.float32,
        **kwargs,
    ) -> ArrayLike:
        """Returns an ArrayLike object copied or referenced from the current
           compute device.

        ---------
        Arguments
        ---------
        buffer_index (int):
            The buffer index specified in the device kernel.

        dtype (type):
            The output type of the returned ArrayLike object. If dtype differs
            from the kernel's buffer type, conversion and truncation rules
            based on the current platform apply.

        """
        output_contents = self._kernel_arg_buffers[buffer_index]
        if isinstance(output_contents, Metal.AGXBuffer):
            output_contents = output_contents.contents().as_buffer(
                output_contents.length()
            )
            output_contents = np.frombuffer(
                output_contents, dtype=dtype
            ).tolist()
        return output_contents
