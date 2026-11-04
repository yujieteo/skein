from __future__ import annotations

import unittest

from replay_diff import differences


def results(value: str) -> str:
    return (
        "verify[2]{stage,id,result,value,limit}:\n"
        '  post,equilibrium,pass,"2.0e-6 N","2.0e-3 N"\n'
        f'  post,expect,pass,"{value}","300 MPa within 5%"\n'
    )


class ExpectTolerance(unittest.TestCase):
    def test_a_change_within_the_tolerance_passes(self) -> None:
        self.assertEqual(differences(results("kirsch: 301.2 MPa"), results("kirsch: 301.3 MPa")), [])

    def test_a_change_beyond_the_tolerance_fails(self) -> None:
        self.assertEqual(len(differences(results("kirsch: 301.2 MPa"), results("kirsch: 330 MPa"))), 1)

    def test_a_changed_name_fails(self) -> None:
        self.assertEqual(len(differences(results("kirsch: 301.2 MPa"), results("other: 301.2 MPa"))), 1)


if __name__ == "__main__":
    unittest.main()
