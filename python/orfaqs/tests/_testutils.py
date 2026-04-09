import pytest
from orfaqs.lib.utils.computeutils import ComputeUtils

accelerators = pytest.mark.skipif(
    not ComputeUtils.compute_accelerator_available(),
    reason='No supported compute accelerator devices were found.',
)
