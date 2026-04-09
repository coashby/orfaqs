from orfaqs.lib.compute.compute_accelerator import ComputeAccelerator
from orfaqs.lib.utils.platformutils import PlatformUtils

if PlatformUtils.is_macos():
    from orfaqs.lib.compute.metal_compute_accelerator import (
        MetalComputeAccelerator,
    )


class ComputeUtils:
    """ComputeUtils"""

    @staticmethod
    def get_default_compute_accelerator() -> ComputeAccelerator | None:
        if (
            PlatformUtils.is_macos()
            and MetalComputeAccelerator.platform_supported()
        ):
            return MetalComputeAccelerator()

        return None

    @staticmethod
    def is_metal_compute_accelerator(
        compute_accelerator: ComputeAccelerator,
    ) -> bool:
        if not PlatformUtils.is_macos():
            return False

        return isinstance(compute_accelerator, MetalComputeAccelerator)

    @staticmethod
    def compute_accelerator_available() -> bool:
        return isinstance(
            ComputeUtils.get_default_compute_accelerator(),
            ComputeAccelerator,
        )
