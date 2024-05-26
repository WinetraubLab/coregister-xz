import math
import numpy as np
import numpy.testing as npt
import unittest
from fit_plane import FitPlane
import plot_fit_plane


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
        
        self.fit_default_plane()
                
    def fit_default_plane(self):
        self.fp = FitPlane(
            self.default_fluorescence_image_points_on_line_pix, 
            self.default_photobleach_line_position_mm, 
            self.default_photobleach_line_group,
            )
            
    def assertAlmostEqualRelative(self, first, second, rel_tol=1e-3, msg=None):
        if not math.isclose(first, second, rel_tol=rel_tol):
            standard_msg = f'{first} != {second} within {rel_tol} relative tolerance'
            self.fail(self._formatMessage(msg, standard_msg))
                
    def test_main_function_runs(self):
        self.fit_default_plane()
    
    def test_example_pixel_size_u(self):
        self.assertAlmostEqual(self.fp.u_norm_mm()*1000, 1, places=0)
   
    def test_example_pixel_size_v(self):
        self.assertAlmostEqual(self.fp.v_norm_mm()*1000, 1, places=0)
        
    def test_example_u_vector_direction(self):
        u_direction = self.fp.u_direction()
        xy_angle = np.degrees(np.arccos(np.dot(u_direction,np.array([1, 0, 0]))))
        z_angle = np.degrees(np.arccos(np.dot(u_direction,np.array([0, 0, 1]))))
        
        self.assertAlmostEqual(xy_angle,45, places=0)
        self.assertAlmostEqual(z_angle ,90, places=1)
        
    def test_plane_normal_computed_correctly(self):
        f = FitPlane([1,0,0],[0,0,1],[10,0,0],method='u,v,h directly')
        n = f.normal_direction()
        
        error_angle = np.degrees(np.arccos(np.dot(n,np.array([0, -1, 0]))))
        self.assertAlmostEqual(error_angle, 0, places=1)

    def test_plane_equation(self):
        # Test with plane z=50
        f = FitPlane([1,0,0],[0,1,0],[10,0,50],method='u,v,h directly')
        a,b,c,d = f.plane_equation()
        
        self.assertAlmostEqual(a, 0, places=1)
        self.assertAlmostEqual(b, 0, places=1)
        self.assertAlmostEqual(c, 1, places=1)
        self.assertAlmostEqual(d, -50, places=1)
    
    def test_example_v_vector_direction(self):
        v_direction = self.fp.v_direction()
        z_angle = np.degrees(np.arccos(np.dot(v_direction,np.array([0, 0, 1]))))
        
        self.assertAlmostEqual(z_angle, 0, places=1)

    def test_u_v_orthogonal_check(self):
        with self.assertRaises(ValueError) as context:
            f = FitPlane([1,0,0],[1,0,0],[10,0,0],method='u,v,h directly')
    
    def test_u_v_nrom_check(self):
        with self.assertRaises(ValueError) as context:
            f = FitPlane([1,0,0],[0,0,2],[10,0,0],method='u,v,h directly')
    
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
                self.default_photobleach_line_group,
                method='points on photobleach lines')
    
    def test_conversion_pixel_position_to_physical_h(self):
        pos_xyz = self.fp.get_xyz_from_uv([0,0])
        self.assertAlmostEqual(pos_xyz[0], self.fp.h[0], places=1)
        self.assertAlmostEqual(pos_xyz[1], self.fp.h[1], places=1)
        self.assertAlmostEqual(pos_xyz[2], self.fp.h[2], places=1)
    
    def test_conversion_pixel_position_to_physical_u(self):
        pos_xyz_1 = self.fp.get_xyz_from_uv([0,0])
        pos_xyz_2 = self.fp.get_xyz_from_uv([1,0])
        self.assertAlmostEqual(pos_xyz_2[0] - pos_xyz_1[0], self.fp.u[0], places=1)
        self.assertAlmostEqual(pos_xyz_2[1] - pos_xyz_1[1], self.fp.u[1], places=1)
        self.assertAlmostEqual(pos_xyz_2[2] - pos_xyz_1[2], self.fp.u[2], places=1)
        
    def test_conversion_pixel_position_to_physical_v(self):
        pos_xyz_1 = self.fp.get_xyz_from_uv([0,0])
        pos_xyz_2 = self.fp.get_xyz_from_uv([0,1])
        self.assertAlmostEqual(pos_xyz_2[0] - pos_xyz_1[0], self.fp.v[0], places=1)
        self.assertAlmostEqual(pos_xyz_2[1] - pos_xyz_1[1], self.fp.v[1], places=1)
        self.assertAlmostEqual(pos_xyz_2[2] - pos_xyz_1[2], self.fp.v[2], places=1)
    
    def test_conversion_pixel_position_to_physical_and_back(self):
        u = 12
        v = 20
        pos_xyz = self.fp.get_xyz_from_uv([u,v])
        pos_uv = self.fp.get_uv_from_xyz(pos_xyz)
        pos_xyz2 = self.fp.get_xyz_from_uv(pos_uv)
        
        self.assertAlmostEqual(pos_uv[0], u, places=1)
        self.assertAlmostEqual(pos_uv[1], v, places=1)
        
        self.assertAlmostEqualRelative(pos_xyz[0], pos_xyz2[0] , rel_tol=0.1) # Accurate up to 10%
        self.assertAlmostEqualRelative(pos_xyz[1], pos_xyz2[1] , rel_tol=0.1) # Accurate up to 10%
        self.assertAlmostEqualRelative(pos_xyz[2], pos_xyz2[2] , rel_tol=0.1) # Accurate up to 10%
        
    def test_check_000_and_physical_conversion(self):
        pos_uv = self.fp.get_uv_from_xyz([0,0,0])
        pos_xyz = self.fp.get_xyz_from_uv(pos_uv)
        pos_xyz = pos_xyz/np.linalg.norm(pos_xyz)
        
        norm_angle = np.degrees(np.arccos(np.dot(pos_xyz, -self.fp.normal_direction())))
        self.assertAlmostEqual(norm_angle, 0, places=1)
        
    def test_conversion_physical_to_pixel_when_point_is_off_plane(self):
        u = 12
        v = 20
        pos_xyz = self.fp.get_xyz_from_uv([u,v])
        pos_xyz = pos_xyz + self.fp.normal_direction()*2 # Move point away from plane
        pos_uv = self.fp.get_uv_from_xyz(pos_xyz)
        
        self.assertAlmostEqual(pos_uv[0], u, places=1)
        self.assertAlmostEqual(pos_uv[1], v, places=1)
        
    def test_recomended_center_pix(self):
        # Modify the pen & paper solution such that lines are centered around origin.
        # In this case, recommended center should be just above origin (0,0,0)
        self.default_photobleach_line_position_mm = [ -0.1, 0, 0.1, 0.2, 0.1, 0, -0.2]
        self.fit_default_plane()
    
        u_coordinates = np.array([item[0] for sublist in self.default_fluorescence_image_points_on_line_pix for item in sublist])
        v_coordinates = np.array([item[1] for sublist in self.default_fluorescence_image_points_on_line_pix for item in sublist])
        
        self.assertAlmostEqualRelative(self.fp.recomended_center_pix[0], np.mean(u_coordinates) , rel_tol=0.1) # Accurate up to 10%
        self.assertAlmostEqualRelative(self.fp.recomended_center_pix[1], np.mean(v_coordinates) , rel_tol=0.1) # Accurate up to 10%

    def verify_two_2D_vectors_point_in_same_direction(self,v1,v2):
        
        # Gather 2D version of each vector
        v1_hat = v1[:2]/np.linalg.norm(v1[:2])
        v2_hat = v2[:2]/np.linalg.norm(v2[:2])
        
        # Compute angle between and make sure it's almost 0
        angle = np.degrees(np.arccos(np.dot(v1_hat,v2_hat)))
        self.assertAlmostEqual(angle, 0, places=0)
        
    
    def test_xy_projection(self):
        pt1, pt2 = self.fp.get_xy_projection()
        
        self.assertEqual(len(pt1),2) # make sure returning values are only x,y
        self.assertEqual(len(pt2),2) # make sure returning values are only x,y
        
        # Verify that pt1-->pt2 is at the same direction as u vector
        self.verify_two_2D_vectors_point_in_same_direction(pt2-pt1, self.fp.u)
    
    def test_xy_projection_with_x_limits(self):
        x_min_mm = -1
        x_max_mm = 1
        
        pt1, pt2 = self.fp.get_xy_projection(min_x_mm=x_min_mm, max_x_mm=x_max_mm)
        pt_min_x, pt_max_x = min(pt1[0], pt2[0]), max(pt1[0], pt2[0])
        
        # Check that x limits are "respected"
        self.assertAlmostEqualRelative(x_min_mm, pt_min_x, rel_tol=0.01)
        self.assertAlmostEqualRelative(x_max_mm, pt_max_x, rel_tol=0.01)
        
        # Verify that pt1-->pt2 is at the same direction as u vector
        self.verify_two_2D_vectors_point_in_same_direction(pt2-pt1, self.fp.u)
        
    def test_xy_projection_with_y_limits(self):
        y_min_mm = -1
        y_max_mm = 1
        
        pt1, pt2 = self.fp.get_xy_projection(min_y_mm=y_min_mm, max_y_mm=y_max_mm)
        pt_min_y, pt_max_y = min(pt1[1], pt2[1]), max(pt1[1], pt2[1])
        
        # Check that x limits are "respected"
        self.assertAlmostEqualRelative(y_min_mm, pt_min_y, rel_tol=0.01)
        self.assertAlmostEqualRelative(y_max_mm, pt_max_y, rel_tol=0.01)
        
        # Verify that pt1-->pt2 is at the same direction as u vector
        self.verify_two_2D_vectors_point_in_same_direction(pt2-pt1, self.fp.u)
        
    def test_plot_fit_plane_doesnt_throw_an_error(self):
        plot_fit_plane.plot_fit_plane(self.fp,[ 0.6, 0.7, 0.8],[-0.2, -0.3, -0.4, -0.6])
        
if __name__ == '__main__':
    unittest.main()
