from orfaqs.lib.compute.compute_accelerator import ComputeAccelerator
from orfaqs.lib.compute.metal_compute_accelerator import (
    MetalComputeAccelerator,
)


class ComputeUtils:
    """ComputeUtils"""

    @staticmethod
    def get_default_compute_accelerator() -> ComputeAccelerator | None:
        if MetalComputeAccelerator.platform_supported():
            return MetalComputeAccelerator()

        return None

    @staticmethod
    def is_metal_compute_accelerator(
        compute_accelerator: ComputeAccelerator,
    ) -> bool:
        return isinstance(compute_accelerator, MetalComputeAccelerator)
