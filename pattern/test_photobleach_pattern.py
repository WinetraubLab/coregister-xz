import math
import numpy as np
import numpy.testing as npt
import matplotlib.pyplot as plt
import unittest
from pattern.photobleach_pattern import Pattern


class TestPhotobleachPattern(unittest.TestCase):

    #def setUp(self):
       
    def test_pattern_not_sorted(self):
        with self.assertRaises(ValueError) as context:
            p = Pattern(v_line_positions=[1,3,2])
        with self.assertRaises(ValueError) as context:
            p = Pattern(v_line_positions=[3,2,1])
        with self.assertRaises(ValueError) as context:
            p = Pattern(h_line_positions=[3,2,1])

    def test_get_number_of_lines(self):
        p = Pattern(v_line_positions=[1,2], h_line_positions=[1,2,3])
        assert(p.get_number_of_lines() == 5)

    def test_project_pattern_onto_plane(self):
        # Test that all projection options work and they give different results
        p = Pattern()
        photobleach_line_position_mm1, photobleach_line_group1 = p.project_pattern_onto_plane(1)
        photobleach_line_position_mm2, photobleach_line_group2 = p.project_pattern_onto_plane(2)
        photobleach_line_position_mm3, photobleach_line_group3 = p.project_pattern_onto_plane(3)
        photobleach_line_position_mm4, photobleach_line_group4 = p.project_pattern_onto_plane(4)

        # Verify result length is the same
        assert(len(photobleach_line_group1) == len(photobleach_line_position_mm1))
        assert(len(photobleach_line_group1) == len(photobleach_line_group2))
        assert(len(photobleach_line_group1) == len(photobleach_line_group3))
        assert(len(photobleach_line_group1) == len(photobleach_line_group4))
        assert(len(photobleach_line_position_mm1) == len(photobleach_line_position_mm2))
        assert(len(photobleach_line_position_mm1) == len(photobleach_line_position_mm3))
        assert(len(photobleach_line_position_mm1) == len(photobleach_line_position_mm3))
        assert(len(photobleach_line_position_mm1) == len(photobleach_line_position_mm4))
        
        # Verify all results are different
        assert(
            (photobleach_line_group1 != photobleach_line_group2) or
            np.any(photobleach_line_position_mm1 != photobleach_line_position_mm2)
            )
        assert(
            (photobleach_line_group1 != photobleach_line_group3) or
            np.any(photobleach_line_position_mm1 != photobleach_line_position_mm3)
            )
        assert(
            (photobleach_line_group1 != photobleach_line_group4) or
            np.any(photobleach_line_position_mm1 != photobleach_line_position_mm4)
            )
        assert(
            (photobleach_line_group2 != photobleach_line_group3) or
            np.any(photobleach_line_position_mm2 != photobleach_line_position_mm3)
            )
        assert(
            (photobleach_line_group2 != photobleach_line_group4) or
            np.any(photobleach_line_position_mm2 != photobleach_line_position_mm4)
            )
        assert(
            (photobleach_line_group3 != photobleach_line_group4) or
            np.any(photobleach_line_position_mm3 != photobleach_line_position_mm4)
            )
        
            
   
if __name__ == '__main__':
    unittest.main()
