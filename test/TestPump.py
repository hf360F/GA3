import unittest

import numpy as np

from main import K as K


class TestPump(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        from main import Pump
        cls.cold = Pump(Pump.COLD)
        cls.hot = Pump(Pump.HOT)

    def test_creation(self):
        self.assertIsNotNone(self.cold)
        self.assertIsNotNone(self.hot)

    def test_float(self):
        self.assertIsNot(self.cold.dp(0.0001), np.nan)
        self.assertIsNot(self.hot.dp(0.0001), np.nan)

    def dont_test_array(self):
        self.assertIsNot(self.cold.dp(np.array([0.0001, 0.0002])), np.nan)
        self.assertIsNot(self.hot.dp(np.array([0.0001, 0.0002])), np.nan)


class TestK2D(unittest.TestCase):
    def test_float(self):
        self.assertIsNot(K(0, 0), np.nan)
        self.assertIsNot(K(0, 25000), np.nan)
        self.assertIsNot(K(0.99, 0), np.nan)
        self.assertIsNot(K(0.99, 25000), np.nan)
        self.assertIsNot(K(0.114, 25.07), np.nan)

    def test_array(self):
        self.assertNotIn(np.nan,K([0,0.5],[0,25000]))


if __name__ == '__main__':
    unittest.main()
