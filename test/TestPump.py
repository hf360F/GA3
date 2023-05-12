import unittest

import numpy as np


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
        self.assertEqual(self.cold.dp(0.0001), 61664.92165351563)
        self.assertEqual(self.hot.dp(0.0001), 48664.711765035114)

    def test_array(self):
        self.assertIsNotNone(self.cold.dp(np.array([0.0001, 0.0002])))
        self.assertIsNotNone(self.hot.dp(np.array([0.0001, 0.0002])))


if __name__ == '__main__':
    unittest.main()
