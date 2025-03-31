'''
PerfUtils
Contains functions for CPU and GPU monitoring and runtime profiling.
'''
import asitop
import asitop.utils
import GPUtil
import logging
import os
import pathlib
import platform
import psutil
import subprocess
import time

from abc import ABC, abstractmethod
from datetime import datetime
from enum import Enum
from threading import (
    Lock,
    Thread
)

from orfaqs.lib.utils.directoryutils import DirectoryUtils

_logger = logging.getLogger(__name__)


class MetricSample:
    def __init__(self,
                 timestamp: int,
                 value: any,
                 value_name: str = None):
        self._timestamp = timestamp
        self._value = value
        self._value_name = 'value'
        if (isinstance(value_name, str) and
                (len(value_name.replace(' ')) > 0)):
            self._value_name = value_name

    def __str__(self):
        return (f'timestamp: {self._timestamp}, '
                f'{self._value_name}: {self._value}')

    @property
    def timestamp(self):
        return self._timestamp

    @property
    def value(self):
        return self._value

    @property
    def map(self) -> dict:
        return {
            'timestamp': self._timestamp,
            self._value_name: self._value
        }


class ElapsedTimeSample:
    '''ElapsedTimeSample'''

    def __init__(self,
                 start_time: int = None,
                 end_time: int = None):
        self._start_time = start_time
        self._end_time = end_time

    @property
    def start_time(self) -> int:
        return self._start_time

    @property
    def end_time(self) -> int:
        return self._end_time

    @property
    def elapsed_time(self) -> int:
        if ((self._start_time is None) or (self._end_time is None)):
            return None

        return self._end_time - self._start_time

    @start_time.setter
    def start_time(self, value: int):
        if ((self._start_time is not None) or
                (value is None)):
            return
        self._start_time = value

    @end_time.setter
    def end_time(self, value: int):
        if ((self._end_time is not None) or
                (value is None)):
            return

        self._end_time = value


class PerfMonitorDeviceType(Enum):
    '''PerfMonitorDeviceType'''

    CPU_DEVICE = 'cpu'
    GPU_DEVICE = 'gpu'


class PerfProfiler:
    '''PerfProfiler'''

    _DEFAULT_CSV_LOG_OUTPUT_DIRECTORY_PATH = '.perf-profiler/'
    _MAX_NUMBER_RECORDS_PER_LOG_FILE = 10 * 1E3
    _MAX_NUMBER_RECORDS_IN_MEMORY = 100 * 1E3

    def __init__(self, profiler_name: str):
        self._name = profiler_name
        self._elapsed_time_samples: dict[str, list[ElapsedTimeSample]] = {}
        self._number_log_file_writes: int = 0
        self._log_file_init()

    @property
    def log_file_index(self) -> int:
        return int(self._number_log_file_writes /
                   self._MAX_NUMBER_RECORDS_PER_LOG_FILE)

    @property
    def csv_log_output_file_path(self) -> pathlib.Path:
        output_directory = DirectoryUtils.make_path_object(
            self._DEFAULT_CSV_LOG_OUTPUT_DIRECTORY_PATH
        )
        output_directory = output_directory.joinpath(self._name)
        file_name = (f'{self._name}'
                     f'-{self.log_file_index}')

        return output_directory.joinpath(file_name).with_suffix('.csv')

    @property
    def _csv_log_file_header(self):
        return ('function_name, '
                'start_timestamp (ns) ,'
                'end_timestamp (ns), '
                'elapsed_time (ns)')

    def _log_file_init(self):
        DirectoryUtils.mkdir_path(self.csv_log_output_file_path.parent)
        # Previous logs are overwritten.
        with open(self.csv_log_output_file_path,
                  'w', encoding='utf-8') as o_file:
            o_file.write(f'{self._csv_log_file_header}\n')

    def _update_log_file(self,
                         function_name: str,
                         elapsed_time_sample: ElapsedTimeSample):
        with open(self.csv_log_output_file_path,
                  'a', encoding='utf-8') as o_file:
            o_file.write(
                (f'{function_name}, '
                 f'{elapsed_time_sample.start_time}, '
                 f'{elapsed_time_sample.end_time}, '
                 f'{elapsed_time_sample.elapsed_time}'
                 '\n')
            )
            self._number_log_file_writes += 1

    @property
    def elapsed_time_samples(self) -> dict[str, list[ElapsedTimeSample]]:
        return self._elapsed_time_samples

    def function_elapsed_time_samples(
            self,
            key: str) -> list[ElapsedTimeSample]:
        return self._elapsed_time_samples.get(key)

    def _validate_stop_perf_timer_call(self,
                                       function_name: str,
                                       samples_list: list):
        calling_stop_before_start_error_message = (
            f'[ERROR] Calling {self.__class__.__qualname__} '
            f'with function_name: "{function_name}, '
            f'before call to {self.__class__}.start_perf_timer.')
        if (function_name not in samples_list):
            _logger.error(calling_stop_before_start_error_message)
            raise RuntimeError(calling_stop_before_start_error_message)

        elif len(samples_list[function_name]) == 0:
            message = ('[ERROR] An unexpected situation occurred, '
                       'elapsed_time_samples array is empty but '
                       f'function_name: {function_name} exists as a key.')
            _logger.error(message)
            raise RuntimeError(message)

        elapsed_time_sample = samples_list[function_name][-1]
        if elapsed_time_sample.start_time is None:
            message = ('[ERROR] An unexpected situation occurred, '
                       'the last element in elapsed_time_samples has no '
                       f'start time for function_name: {function_name}.')
            _logger.error(message)
            raise RuntimeError(message)
        elif elapsed_time_sample.end_time is not None:
            _logger.error(calling_stop_before_start_error_message)
            raise RuntimeError(calling_stop_before_start_error_message)

    def start_perf_timer(self, function_name: str):
        elapsed_time_sample = ElapsedTimeSample(time.perf_counter_ns())
        if function_name not in self._elapsed_time_samples:
            self._elapsed_time_samples[function_name] = []

        self._elapsed_time_samples[function_name].append(
            elapsed_time_sample
        )

    def stop_perf_timer(self, function_name: str):
        self._validate_stop_perf_timer_call(
            function_name,
            self._elapsed_time_samples
        )
        self._elapsed_time_samples[function_name][-1].end_time = (
            time.perf_counter_ns()
        )
        self._update_log_file(
            function_name,
            self._elapsed_time_samples[function_name][-1]
        )


class _PerfMonitor(ABC):
    '''_PerfMonitor'''

    _DEFAULT_POLLING_RATE = 1.0
    _DEFAULT_CSV_LOG_OUTPUT_DIRECTORY_PATH = '.perf-monitor/'
    _MAX_NUMBER_RECORDS_PER_LOG_FILE = 10 * 1E3
    _MAX_NUMBER_RECORDS_IN_MEMORY = 100 * 1E3
    _VALUE_PRECISION_DECIMAL_LIMIT = 2

    def __del__(self):
        self.stop_metrics_polling()

    def __init__(self, monitor_name: str, device_type: PerfMonitorDeviceType):
        self._name = monitor_name
        self._device_type = device_type
        self._device_name: str = None
        self._device_type = device_type
        self._number_cores: int = None
        self._compute_utilization_samples: list[MetricSample] = []
        self._allocated_memory_bytes_samples: list[MetricSample] = []
        self._total_memory: int = None
        self._metrics_polling_thread: Thread = None
        self._metrics_data_lock: Lock = Lock()
        self._poll_metrics: bool = False
        self._polling_rate: float = _PerfMonitor._DEFAULT_POLLING_RATE
        self._creation_datetime_str = datetime.now().strftime('%Y%m%d-%H%M%S')
        self._number_log_file_writes: int = 0
        self._log_file_init()

    @property
    def log_file_index(self):
        return int(self._number_log_file_writes /
                   self._MAX_NUMBER_RECORDS_PER_LOG_FILE)

    @property
    def csv_log_output_file_path(self) -> pathlib.Path:
        output_directory = DirectoryUtils.make_path_object(
            self._DEFAULT_CSV_LOG_OUTPUT_DIRECTORY_PATH
        )
        output_directory = output_directory.joinpath(
            self._creation_datetime_str
        )
        class_name = self.__class__.__name__.lower()
        file_name = (f'{class_name}'
                     f'-{self._name}'
                     f'-{self._creation_datetime_str}'
                     f'-{self.log_file_index}')

        return output_directory.joinpath(file_name).with_suffix('.csv')

    @property
    def _csv_log_file_header(self):
        return ('timestamp (ns), '
                f'{self._device_type.value}_utilization, '
                f'{self._device_type.value}_memory_utilization')

    def _log_file_init(self):
        output_file_path = self.csv_log_output_file_path
        DirectoryUtils.mkdir_path(output_file_path.parent)
        with open(output_file_path, 'w', encoding='utf-8') as o_file:
            o_file.write(f'{self._csv_log_file_header}\n')

    def _update_log_file(self,
                         timestamp,
                         compute_utilization,
                         memory_utilization):
        with open(self.csv_log_output_file_path,
                  'a', encoding='utf-8') as o_file:
            o_file.write(
                (f'{timestamp}, '
                 f'{compute_utilization}, '
                 f'{memory_utilization}'
                 '\n')
            )
            self._number_log_file_writes += 1

    @property
    def device_name(self) -> str:
        return self._device_name

    @property
    def number_cores(self) -> int:
        return self._number_cores

    @property
    def compute_utilization_samples(self) -> list[MetricSample]:
        utilization_samples: list[MetricSample] = None
        with self._metrics_data_lock:
            utilization_samples = self._compute_utilization_samples.copy()

        return utilization_samples

    @property
    def allocated_memory_bytes_samples(self) -> list[MetricSample]:
        memory_samples: list[MetricSample] = None
        with self._metrics_data_lock:
            memory_samples = self._allocated_memory_bytes_samples.copy()

        return memory_samples

    @property
    def memory_utilization_samples(self) -> list[MetricSample]:
        utilization_samples: list[MetricSample] = {}
        with self._metrics_data_lock:
            for sample in self._allocated_memory_bytes_samples:
                utilization = round(
                    ((100.0 * sample.value) / self._total_memory),
                    self._VALUE_PRECISION_DECIMAL_LIMIT
                )
                utilization_samples.append(
                    MetricSample(sample.timestamp, utilization)
                )

        return utilization_samples

    @property
    def total_memory(self) -> int:
        return self._total_memory

    def last_compute_utilization_sample(self) -> MetricSample:
        if (isinstance(self._compute_utilization_samples, list) and
                (len(self._compute_utilization_samples) > 0)):
            with self._metrics_data_lock:
                return self._compute_utilization_samples[-1]

        return None

    def last_memory_utilization_sample(self) -> MetricSample:
        if (isinstance(self._allocated_memory_bytes_samples, list) and
                (len(self._allocated_memory_bytes_samples) > 0)):
            with self._metrics_data_lock:
                sample = self._allocated_memory_bytes_samples[-1]
                utilization = round(
                    ((100.0 * sample.value) / self._total_memory),
                    self._VALUE_PRECISION_DECIMAL_LIMIT
                )

                return MetricSample(sample.timestamp, utilization)

        return None

    def _update_compute_utilization_samples(self,
                                            timestamp: int,
                                            value: float):
        if value is not None:
            value = round(value, self._VALUE_PRECISION_DECIMAL_LIMIT)
        with self._metrics_data_lock:
            self._compute_utilization_samples.append(
                MetricSample(timestamp, value)
            )
            if (len(self._compute_utilization_samples) >
                    self._MAX_NUMBER_RECORDS_IN_MEMORY):
                self._compute_utilization_samples = (
                    self._compute_utilization_samples[1:]
                )

    def _update_allocated_memory_bytes_samples(self,
                                               timestamp: int,
                                               value: float):
        if value is not None:
            value = round(value, self._VALUE_PRECISION_DECIMAL_LIMIT)
        with self._metrics_data_lock:
            self._allocated_memory_bytes_samples.append(
                MetricSample(timestamp, value)
            )
            if (len(self._allocated_memory_bytes_samples) >
                    self._MAX_NUMBER_RECORDS_IN_MEMORY):
                self._allocated_memory_bytes_samples = (
                    self._allocated_memory_bytes_samples[1:]
                )

    @abstractmethod
    def _collect_metrics(self) -> tuple[int, float, int]:
        '''
        Returns the processing unit's utilization and memory usage as a dict.
        '''
        pass

    def start_metrics_polling(self, nop_time_seconds: float = None):
        if isinstance(self._metrics_polling_thread, Thread):
            return

        def _metrics_polling():
            self._poll_metrics = True
            while self._poll_metrics:
                (timestamp,
                 compute_utilization,
                 allocated_memory_bytes) = self._collect_metrics()

                self._update_compute_utilization_samples(
                    timestamp,
                    compute_utilization
                )

                self._update_allocated_memory_bytes_samples(
                    timestamp,
                    allocated_memory_bytes
                )

                memory_utilization = None
                if allocated_memory_bytes is not None:
                    memory_utilization = (
                        (100.0 * allocated_memory_bytes) / self._total_memory
                    )
                    memory_utilization = round(
                        memory_utilization,
                        self._VALUE_PRECISION_DECIMAL_LIMIT
                    )

                compute_utilization = round(
                    compute_utilization,
                    self._VALUE_PRECISION_DECIMAL_LIMIT
                )
                self._update_log_file(
                    timestamp,
                    compute_utilization,
                    memory_utilization
                )
                time.sleep(self._polling_rate)

        self._metrics_polling_thread = Thread(target=_metrics_polling)
        self._metrics_polling_thread.start()
        # Optional wait time before returning to the caller.
        # This allows the profiler to establish some baseline measurements.
        if nop_time_seconds is not None:
            time.sleep(nop_time_seconds)

    def stop_metrics_polling(self):
        if self._metrics_polling_thread is None:
            return
        elif isinstance(self._metrics_polling_thread, Thread):
            self._poll_metrics = False
            if self._metrics_polling_thread.is_alive():
                self._metrics_polling_thread.join()
            self._metrics_polling_thread = None

    def export_metrics_as_csv(self, output_path: str | os.PathLike):
        pass

    def export_metrics_as_json(self, output_path: str | os.PathLike):
        pass


class CPUPerfMonitor(_PerfMonitor):
    '''CpuPerfMonitor'''

    def __init__(self, monitor_name):
        super().__init__(
            monitor_name,
            device_type=PerfMonitorDeviceType.CPU_DEVICE
        )

        self._total_memory = psutil.virtual_memory().total

    def _collect_metrics(self) -> tuple[int, float, float]:
        '''
        Returns the processing unit's utilization and memory usage.
        return (timestamp, compute_utilization, allocated_memory_bytes)
        '''
        timestamp = time.perf_counter_ns()
        compute_utilization = psutil.cpu_percent()
        psutil_memory_utilization = psutil.virtual_memory()
        allocated_memory_bytes = None
        if platform.system() == "Windows":
            allocated_memory_bytes = psutil_memory_utilization.used
        else:
            allocated_memory_bytes = psutil_memory_utilization.active
        return (timestamp, compute_utilization, allocated_memory_bytes)

    def collect_metrics(self) -> dict[any]:
        '''
        Returns the processing unit's utilization and memory usage as a dict.
        '''
        timestamp = time.perf_counter_ns()
        psutil_cpu_utilization = psutil.cpu_times_percent()
        psutil_memory_utilization = psutil.virtual_memory()
        allocated_memory_bytes = None
        if platform.system() == "Windows":
            allocated_memory_bytes = psutil_memory_utilization.used
        else:
            allocated_memory_bytes = psutil_memory_utilization.active
        return {
            'timestamp': timestamp,
            'cpu_utilization': {
                'user': psutil_cpu_utilization.user,
                'system': psutil_cpu_utilization.system,
                'idle': psutil_cpu_utilization.idle
            },
            'memory_utilization': {
                'total_memory': psutil_memory_utilization.total,
                'allocated_memory_bytes': allocated_memory_bytes,
                'available_memory': psutil_memory_utilization.available
            }
        }


class AppleSiGpuPerfMonitor(_PerfMonitor):
    '''AppleSiGpuPerfMonitor'''

    def __init__(self, monitor_name, gpu_index: int = 0):
        super().__init__(
            monitor_name,
            device_type=PerfMonitorDeviceType.GPU_DEVICE
        )

        self._has_gpus = asitop.utils.get_gpu_cores() > 0
        self._number_gpus = asitop.utils.get_gpu_cores()
        self._gpu_index = None
        self._default_gpu: GPUtil.GPU = None
        if self._has_gpus:
            if ((gpu_index < 0) or
                    (gpu_index >= self._number_gpus)):
                message = ('[ERROR] Invalid index for GPU selection.\n'
                           f'GPU Index Given: {gpu_index}\n'
                           f'Available GPU Count: {self._number_gpus}')
                _logger.error(message)
                raise ValueError(message)

            self._gpu_index = gpu_index

    @staticmethod
    def gpus_available() -> bool:
        return asitop.utils.get_gpu_cores() > 0

    @property
    def number_gpus(self) -> int:
        return self._number_gpus

    @staticmethod
    def _gpu_utilization() -> float:
        command = [
            'sudo',
            'powermetrics',
            '--samplers',
            'gpu_power',
            '-i',
            '500',
            '-n',
            '1'
        ]
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        (stdout, _) = process.communicate()
        metrics = stdout.decode('utf-8')
        metrics_lines = metrics.splitlines()
        utilization = None
        for line in metrics_lines:
            key = 'GPU HW active residency:'
            if key in line:
                line = line.replace(key, '')
                utilization_str = line.strip().split(' ')[0]
                utilization = float(utilization_str.replace('%', '').strip())
                break

        return utilization

    def _collect_metrics(self) -> tuple[int, float, float]:
        '''
        Returns the processing unit's utilization and memory usage.
        return (timestamp, compute_utilization, allocated_memory_bytes)
        '''
        timestamp = time.perf_counter_ns()
        compute_utilization = self._gpu_utilization()
        allocated_memory_bytes = None
        return (timestamp, compute_utilization, allocated_memory_bytes)

    def collect_metrics(self) -> dict[any]:
        '''
        Returns the processing unit's utilization and memory usage as a dict.
        '''
        timestamp = time.perf_counter_ns()
        gpu_utilization = self._gpu_utilization()
        return {
            'timestamp': timestamp,
            'gpu_utilization': {
                'active': gpu_utilization,
                'idle': (100 - gpu_utilization)
            },
            'memory_utilization': {
            }
        }


class NvidiaGpuPerfMonitor(_PerfMonitor):
    '''NvidiaGpuPerfMonitor'''

    def __init__(self, monitor_name, gpu_index: int = 0):
        super().__init__(
            monitor_name,
            device_type=PerfMonitorDeviceType.GPU_DEVICE
        )

        self._available_gpus = GPUtil.getGPUs()
        self._number_gpus = len(self._available_gpus)
        self._has_gpus = self._number_gpus > 0
        self._gpu_index = None
        self._default_gpu: GPUtil.GPU = None
        if self._has_gpus:
            if ((gpu_index < 0) or
                    (gpu_index >= self._number_gpus)):
                message = ('[ERROR] Invalid index for GPU selection.\n'
                           f'GPU Index Given: {gpu_index}\n'
                           f'Available GPU Count: {self._number_gpus}')
                _logger.error(message)
                raise ValueError(message)

            self._gpu_index = gpu_index
            self._default_gpu = self._available_gpus[self._gpu_index]
            self._total_memory = self._default_gpu.memoryTotal

    @staticmethod
    def gpus_available() -> bool:
        return len(GPUtil.getGPUs()) > 0

    @property
    def number_gpus(self) -> int:
        return self._number_gpus

    @property
    def default_gpu(self) -> GPUtil.GPU:
        if not self.gpus_available():
            return None
        return GPUtil.getGPUs()[0]

    def _collect_metrics(self) -> tuple[int, float, float]:
        '''
        Returns the processing unit's utilization and memory usage.
        return (timestamp, compute_utilization, allocated_memory_bytes)
        '''
        if self.default_gpu is None:
            return
        timestamp = time.perf_counter_ns()
        compute_utilization = self.default_gpu.load * 100
        allocated_memory_bytes = self.default_gpu.memoryUsed * 1024 * 1024
        return (timestamp, compute_utilization, allocated_memory_bytes)

    def collect_metrics(self) -> dict[any]:
        '''
        Returns the processing unit's utilization and memory usage as a dict.
        '''
        if self.default_gpu is None:
            return
        timestamp = time.perf_counter_ns()
        compute_utilization = self.default_gpu.load * 100
        total_memory_bytes = self.default_gpu.memoryTotal * 1024 * 1024
        allocated_memory_bytes = self.default_gpu.memoryUsed * 1024 * 1024
        return {
            'timestamp': timestamp,
            'gpu_utilization': compute_utilization,
            'memory_utilization': {
                'total_memory_bytes': total_memory_bytes,
                'allocated_memory_bytes': allocated_memory_bytes
            }
        }


class GPUPerfMonitorSelector:
    '''GpuPerfMonitorSelector'''

    @staticmethod
    def load_platform_gpu(monitor_name, gpu_index: int = 0):
        if NvidiaGpuPerfMonitor.gpus_available():
            return NvidiaGpuPerfMonitor(monitor_name, gpu_index)
        elif AppleSiGpuPerfMonitor.gpus_available():
            return AppleSiGpuPerfMonitor(monitor_name, gpu_index)

        return None
