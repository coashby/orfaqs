"""
Orfaqs App
"""

from orfaqs.lib.utils.directoryutils import DirectoryUtils

from abc import ABC, abstractmethod

from orfaqs.lib.utils.jsonutils import JsonUtils
from orfaqs.lib.utils.perfutils import PerfProfiler


class ORFaqsApp(ABC):
    """OrfaqsApp"""

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
    def launch_json_option_name() -> str:
        return 'launch_json'

    @staticmethod
    def process_ui_kwargs(
        launch_json_key: str = None,
        **kwargs,
    ) -> dict[str, any]:
        if launch_json_key is None:
            launch_json_key = ORFaqsApp.launch_json_option_name()

        if kwargs.get(launch_json_key) is not None:
            launch_args: dict[str, any] = {}
            launch_json_path = kwargs[launch_json_key]
            launch_args = JsonUtils.read_json(launch_json_path)
            # Only replace input values that have not been set.
            for key, value in launch_args.items():
                if kwargs.get(key) is None:
                    kwargs[key] = value
                elif kwargs.get(key) is False:
                    kwargs[key] = value

        return kwargs

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
