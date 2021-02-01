import numpy as np
import unittest
from mg1_model import MG1Queue


class MG1Tests(unittest.TestCase):
    def setUp(self) -> None:
        theta1 = 1
        theta2 = 3
        theta3 = 0.4
        self.model = MG1Queue([theta1, theta2, theta3])
        self.rng = np.random.RandomState(seed=42)

    def test_check_input(self):
        self.assertTrue(not self.model._check_input([1, 2, -1]))
        self.assertTrue(not self.model._check_input([-1, 2, 1]))
        self.assertTrue(not self.model._check_input([3, 2, 1]))
        self.assertTrue(not self.model._check_input([1, 0, 1]))

    def test_forward_sim(self):
        out = self.model.forward_simulate([1, 3, 0.4], num_forward_simulations=2, rng=self.rng)
        self.assertTrue(np.allclose(out[0], np.array([4.07459884, 2.58775259, 1.31198904, 2.73235229, 2.41614516])))
        print(out)
