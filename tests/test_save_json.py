import json

import numpy as np

from res_reliability.analyses import save_json


def test_save_json_writes_null_for_non_finite_numbers(tmp_path):
    save_json("x.json", {"a": float("nan"), "b": [np.float64(np.inf), 1.5], "c": np.float32(2.0)}, tmp_path)
    text = (tmp_path / "x.json").read_text(encoding="utf-8")
    assert "NaN" not in text and "Infinity" not in text
    assert json.loads(text) == {"a": None, "b": [None, 1.5], "c": 2.0}
