'''
Orfaqs App
'''
from orfaqs.lib.utils.directoryutils import DirectoryUtils

from abc import (
    ABC,
    abstractmethod
)

from orfaqs.lib.utils.perfutils import PerfProfiler


class OrfaqsApp(ABC):
    '''OrfaqsApp'''

    @staticmethod
    @abstractmethod
    def program_name() -> str:
        pass

    @classmethod
    def default_output_directory(cls):
        base_path = DirectoryUtils.make_path_object('./.orfaqs-apps/outputs')
        formatted_program_directory = (
            cls.program_name().replace(' ', '-').lower()
        )
        return base_path.joinpath(formatted_program_directory)

    @staticmethod
    @abstractmethod
    def cli():
        pass

    @staticmethod
    @abstractmethod
    def _run(**kwargs):
        pass

    @classmethod
    def run(cls, **kwargs):
        perf_profiler = PerfProfiler(cls.__name__)
        function_name = cls.__class__.__qualname__
        perf_profiler.start_perf_timer(function_name)
        cls._run(**kwargs)
        perf_profiler.stop_perf_timer(function_name)
