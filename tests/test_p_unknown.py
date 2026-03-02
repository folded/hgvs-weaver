import pytest

import weaver


@pytest.mark.parametrize("s", ["NP_037359.3:p.?", "NP_037359.3:p.(=)"])
def test_parse_p_unknown(s: str) -> None:
    v = weaver.parse(s)
    assert v is not None

    d = v.to_dict()
    assert isinstance(d, dict)
    assert str(v) == s
