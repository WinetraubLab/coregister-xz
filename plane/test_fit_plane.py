import numpy as np
import numpy.testing as npt
import unittest
from fit_plane import FitPlane

class TestFitPlane(unittest.TestCase):

    # Define default pattern
    def setUp(self):
        # For this test, we hand computed a few points and make sure the plane fits correctly
        
        # Define the pattern from above
        self.default_photobleach_line_position_mm = [ 0.6, 0.7, 0.8, -0.2, -0.3, -0.4, -0.6] # Offset of each line from the center
        self.default_photobleach_line_group = ['v', 'v', 'v', 'h', 'h', 'h', 'h'] # v means the line is vertical line, h means horizontal line, 's' is surface line
        
        # Identify positions on the fluorescence image that are on each lines
        self.default_fluorescence_image_points_on_line_pix = [
                [ [ 143, 0], [ 143, 1] ], # A few points on line 0 (v)
                [ [ 272, 0], [ 272, 1] ], # A few points on line 1 (v)
                [ [ 412, 0], [ 412, 1] ], # ...
                [ [1359, 0], [1359, 1] ],
                [ [1492, 0], [1492, 1] ],
                [ [1625, 0], [1625, 1] ],
                [ [1894, 0], [1894, 1] ], # A few points on line 6 (h)
                ]
                
    def test_main_function_runs(self):
        self.fp = FitPlane(
            self.default_fluorescence_image_points_on_line_pix, 
            self.default_photobleach_line_position_mm, 
            self.default_photobleach_line_group)
    
    def test_example_pixel_size_u(self):
        self.test_main_function_runs()
        self.assertAlmostEqual(self.fp.u_norm_mm()*1000, 1, places=0)
   
    def test_example_pixel_size_v(self):
        self.test_main_function_runs()
        self.assertAlmostEqual(self.fp.v_norm_mm()*1000, 1, places=0)
        
    def test_example_u_vector_direction(self):
        self.test_main_function_runs()
        u_direction = self.fp.u_direction()
        xy_angle = np.degrees(np.arccos(np.dot(u_direction,np.array([1, 0, 0]))))
        z_angle = np.degrees(np.arccos(np.dot(u_direction,np.array([0, 0, 1]))))
        
        self.assertAlmostEqual(xy_angle,45, places=0)
        self.assertAlmostEqual(z_angle ,90, places=1)
    
    def test_example_v_vector_direction(self):
        self.test_main_function_runs()
        v_direction = self.fp.v_direction()
        z_angle = np.degrees(np.arccos(np.dot(v_direction,np.array([0, 0, 1]))))
        
        self.assertAlmostEqual(z_angle ,0, places=1)    
    
    def test_tilted_image(self):
        # For this test, we will look at the case where fluorescence image is tilted.
        # Our code should return an error
        
        # Define the rotation matrix for 10 degrees
        theta = np.radians(10)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        rotation_matrix = np.array([
            [cos_theta, -sin_theta],
            [sin_theta, cos_theta]
        ])
        
        # Apply the rotation matrix to each point
        rotated_fluorescence_image_points_on_line_pix = self.default_fluorescence_image_points_on_line_pix @ rotation_matrix.T
        
        # Fit a plane
        with self.assertRaises(ValueError) as context:
            fp = FitPlane(
                rotated_fluorescence_image_points_on_line_pix, 
                self.default_photobleach_line_position_mm, 
                self.default_photobleach_line_group)

if __name__ == '__main__':
    unittest.main()
