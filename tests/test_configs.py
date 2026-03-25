from pathlib import Path

import pytest
import quark
from .conftest import TEST_DIR

configs = list((TEST_DIR / "configs").glob("*.yml"))


@pytest.mark.parametrize(
    "config",
    configs,
    ids=[str(config.stem) for config in configs],
)
def test_configs(config: Path) -> None:
    quark.start(["-c", str(config)])
