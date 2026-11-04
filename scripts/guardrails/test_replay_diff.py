from __future__ import annotations

import unittest

from pathlib import Path
from tempfile import TemporaryDirectory

from replay_diff import compare, differences


def results(value: str, equilibrium: str = "2.0e-6 N") -> str:
    return (
        "verify[2]{stage,id,result,value,limit}:\n"
        f'  post,equilibrium,pass,"{equilibrium}","2.0e-3 N"\n'
        f'  post,expect,pass,"{value}","300 MPa within 5%"\n'
    )


class ExpectTolerance(unittest.TestCase):
    def test_a_change_within_the_tolerance_passes(self) -> None:
        self.assertEqual(differences(results("kirsch: 301.2 MPa"), results("kirsch: 301.3 MPa")), [])

    def test_a_change_beyond_the_tolerance_fails(self) -> None:
        self.assertEqual(len(differences(results("kirsch: 301.2 MPa"), results("kirsch: 330 MPa"))), 1)

    def test_a_changed_name_fails(self) -> None:
        self.assertEqual(len(differences(results("kirsch: 301.2 MPa"), results("other: 301.2 MPa"))), 1)

    def test_a_plain_limit_is_a_pass_threshold_not_a_drift(self) -> None:
        self.assertEqual(
            len(differences(results("kirsch: 301.2 MPa"), results("kirsch: 301.2 MPa", "1.9e-3 N"))), 1
        )

    def test_a_plain_limit_allows_only_rounding(self) -> None:
        self.assertEqual(differences(stress("301.2 MPa"), stress("301.2000001 MPa")), [])

    def test_a_count_must_match_exactly(self) -> None:
        self.assertEqual(len(differences(steps("1000000 steps"), steps("1000001 steps"))), 1)


def stress(value: str) -> str:
    return (
        "verify[1]{stage,id,result,value,limit}:\n"
        f'  post,stress,pass,"{value}","301.2 MPa"\n'
    )


def steps(value: str) -> str:
    return (
        "verify[1]{stage,id,result,value,limit}:\n"
        f'  solve,finite,pass,"{value}","each value finite"\n'
    )


class NewModel(unittest.TestCase):
    def run_compare(self, base_names: list[str], head_names: list[str]) -> list[str]:
        with TemporaryDirectory() as tmp:
            base, head = Path(tmp) / "base", Path(tmp) / "head"
            for directory, names in ((base, base_names), (head, head_names)):
                directory.mkdir()
                for name in names:
                    (directory / name).write_text("x")
            return compare(base, head)[1]

    def test_a_model_that_the_base_cannot_build_is_new(self) -> None:
        self.assertEqual(self.run_compare(["m.fail"], ["m.c"]), [])

    def test_a_model_that_the_branch_cannot_build_is_reported(self) -> None:
        self.assertEqual(len(self.run_compare(["m.c"], ["m.fail"])), 1)


if __name__ == "__main__":
    unittest.main()
