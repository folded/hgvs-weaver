import json

import weaver


def test_parse(s):
    try:
        v = weaver.parse(s)
        d = v.to_dict()
        print(f"Parsed {s}: {v}")
        print(json.dumps(d, indent=2))
        return v
    except Exception as e:
        print(f"Failed to parse {s}: {e}")
        return None


test_parse("NP_037359.3:p.?")
test_parse("NP_037359.3:p.(=)")
