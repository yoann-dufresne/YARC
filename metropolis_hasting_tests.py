import unittest
from metropolis_hasting import MetropolisHasting1D


class MetropolisTest(unittest.TestCase):
    """ Tests for the class MetropolisHasting1D """

    def test_good_distrib(self):
        """Input a good sequence for MH initialization."""
        distrib = [(3,1), (4,1), (6,3), (7,5), (12, 2), (17,1)]
        MetropolisHasting1D(distrib, 2, burning_steps=10)

    def test_duplicate_value_distrib(self):
        """Input a bad sequence with duplicate x value
        Must raise an exception
        """
        distrib = [(3,1), (4,1), (6,3), (7,5), (12, 2), (17,1), (17,2)]
        with self.assertRaises(Exception):
            MetropolisHasting1D(distrib, 2, burning_steps=10)

    def test_unordered_value_distrib(self):
        """Input a bad sequence with unordered x values
        Must raise an exception
        """
        distrib = [(3,1), (4,1), (6,3), (7,5), (12, 2), (17,1), (4,2)]
        with self.assertRaises(Exception):
            MetropolisHasting1D(distrib, 2, burning_steps=10)

    def test_model_distribution(self):
        """Init a distribution with steps with round values"""
        distrib = [(3,1), (4,1), (6,3), (7,5), (12, 2), (17,1)]
        metro = MetropolisHasting1D(distrib, 2, burning_steps=10)
        self.assertEqual(metro.distribution,
            [(3,1), (5,2), (7,5), (9, 3.8), (11,2.6), (13, 1.8), (15, 1.4), (17,1)])


if __name__ == '__main__':
    unittest.main()
