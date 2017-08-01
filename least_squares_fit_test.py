import unittest
from StringIO import StringIO
from textwrap import dedent

import newick
import numpy as np

from least_squares_fit import post_order, read_distance_matrix

class DistanceMatrixTest(unittest.TestCase):
    def setUp(self):
        newick_str = '(((A:1, B:2)I:3.5, (C:4, D:5)J:6), E:7)K;'
        self.tree = newick.loads(newick_str)[0]
        self.matrix_str = dedent('''\
        A\tB\t3
        A\tC\t14.5
        B\tD\t16.5
        B\tE\t13.5
        C\tD\t9
        C\tE\t18
        D\tE\t19
        A\tD\t15.5
        A\tE\t12.5
        B\tC\t15.5''')

    def test_distance_matrix_happy(self):
        """Tests the "happy path": makes sure the distance matrix looks right in normal conditions."""
        fake_file = StringIO(self.matrix_str)
        po = post_order(self.tree)
        matrix = read_distance_matrix(fake_file, po)
        np.testing.assert_array_equal(matrix, np.matrix([[0,    3,    14.5, 15.5, 12.5],
                                                         [3,    0,    15.5, 16.5, 13.5],
                                                         [14.5, 15.5, 0,       9, 18],
                                                         [15.5, 16.5, 9,       0, 19],
                                                         [12.5, 13.5, 18,     19,  0]]))

    def test_distance_matrix_not_enough_entries(self):
        """The distance matrix reader should complain if we forget an entry."""
        one_line_missing = "\n".join(self.matrix_str.split("\n")[:-1])
        fake_file = StringIO(one_line_missing)
        po = post_order(self.tree)
        with self.assertRaises(RuntimeError):
            read_distance_matrix(fake_file, po)

    def test_distance_matrix_additional_entries(self):
        """Additional entries are OK in the distance matrix if they match earlier entries."""
        matrix_str = self.matrix_str + dedent('''
        B\tA\t3
        E\tC\t18
        ''')
        fake_file = StringIO(matrix_str)
        po = post_order(self.tree)
        matrix = read_distance_matrix(fake_file, po)
        np.testing.assert_array_equal(matrix, np.matrix([[0,    3,    14.5, 15.5, 12.5],
                                                         [3,    0,    15.5, 16.5, 13.5],
                                                         [14.5, 15.5, 0,       9, 18],
                                                         [15.5, 16.5, 9,       0, 19],
                                                         [12.5, 13.5, 18,     19,  0]]))

    def test_distance_matrix_conflicting_entries(self):
        """Conflicting entries in the input distance matrix should be rejected."""
        matrix_str = self.matrix_str + dedent('''
        B\tA\t10
        E\tC\t300
        ''')
        fake_file = StringIO(matrix_str)
        po = post_order(self.tree)
        with self.assertRaises(RuntimeError):
            read_distance_matrix(fake_file, po)

if __name__ == '__main__':
    unittest.main()
