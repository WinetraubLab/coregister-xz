import math
import numpy as np
import numpy.testing as npt
import matplotlib.pyplot as plt
import unittest
import image_processing.fluorescent_image_preprocessing as t
from image_processing.draw_fluorescent_image import draw_fluorescent_image
import cv2
import matplotlib.pyplot as plt

class TestFluorescentImagePreProcessing(unittest.TestCase):

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

    def test_rotate_pt_list_rand_rotate(self):
        pt0 = self.default_fluorescence_image_points_on_line_pix;
        rotated_pt = t._rotate_pt_list(pt0,np.pi/4);

    def test_rotate_pt_list_zero_rotate(self):
        pt0 = self.default_fluorescence_image_points_on_line_pix;
        rotated_pt = t._rotate_pt_list(pt0,0);

        # Validate shape is the same
        assert(len(rotated_pt) == len(pt0))
        assert(len(rotated_pt[0]) == len(pt0[0]))
        assert(len(rotated_pt[0][0]) == len(pt0[0][0]))

        # Validate some of the values are the same
        self.assertAlmostEqual(pt0[0][0][0], rotated_pt[0][0][0], places=1)
        self.assertAlmostEqual(pt0[0][0][1], rotated_pt[0][0][1], places=1)
        self.assertAlmostEqual(pt0[0][1][0], rotated_pt[0][1][0], places=1)
        self.assertAlmostEqual(pt0[3][1][0], rotated_pt[3][1][0], places=1)

    def test_rotate_pt_list_90deg_rotate(self):
        pt0 = self.default_fluorescence_image_points_on_line_pix;
        rotated_pt = t._rotate_pt_list(pt0,np.pi/2);

        # Validate some of the values are the same
        self.assertAlmostEqual(pt0[0][0][0], -rotated_pt[0][0][1], places=1)
        self.assertAlmostEqual(pt0[0][0][1], rotated_pt[0][0][0], places=1)
        self.assertAlmostEqual(pt0[0][1][0], -rotated_pt[0][1][1], places=1)
        self.assertAlmostEqual(pt0[0][1][1], rotated_pt[0][1][0], places=1)
                   

    def test_rotate_points_to_optimize_uv(self):

        def test_one_angle(angle):
            # Rotate points
            rotated_pt = t._rotate_pt_list(
                self.default_fluorescence_image_points_on_line_pix,
                angle);

            # Optimize, see that the result angle is back to the optimum (0 rotation)
            optimized_theta = t._rotate_points_to_optimize_uv(
                rotated_pt,
                self.default_photobleach_line_group)

            # Check that optimal rotation brings back to default
            self.assertAlmostEqual(optimized_theta,-angle, places=2)

        test_one_angle( 45*np.pi/180)
        test_one_angle(-45*np.pi/180)
        test_one_angle( 10*np.pi/180)

    def test_image_rotation_to_optimize_uv(self):

        pt0 = self.default_fluorescence_image_points_on_line_pix;

        # Rotate points
        rotated_pt = t._rotate_pt_list(
                self.default_fluorescence_image_points_on_line_pix,
                -np.pi/4,
                center_x=50-1,
                center_y=100-1,
                );

        # Set black input image with white pixel
        image = np.zeros((100, 200, 3), dtype=np.uint8)
        image[int(rotated_pt[0][1][1]),int(rotated_pt[0][1][0]),2] = 255;

        # Make sure the pixel does show color
        self.assertGreater(image[int(rotated_pt[0][1][1]),int(rotated_pt[0][1][0]),2], 50);

        # Find optimal solution and rotate back    
        pt_out, image_out = t.rotate_image_to_meet_consistency_assumptions(
            rotated_pt,
            self.default_photobleach_line_group,
            image)

        # Verify image shape didn't change
        assert(image_out.shape[0] == image_out.shape[0])
        assert(image_out.shape[1] == image_out.shape[1])
        assert(image_out.shape[2] == image_out.shape[2])

        if False: # Change to true if you would like to plot
            # Plot original image
            plt.figure(figsize=(10, 5))
            plt.subplot(1, 2, 1)
            plt.title("Original Image")
            plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))  # Convert BGR to RGB for proper color display
            plt.axis('off')

            # Plot rotated image
            plt.subplot(1, 2, 2)
            plt.title("Rotated Image")
            plt.imshow(cv2.cvtColor(image_out, cv2.COLOR_BGR2RGB))  # Convert BGR to RGB for proper color display
            plt.axis('off')
            plt.show()

        # Verify that the rotated pixel is white
        self.assertGreater(image_out[int(pt_out[0][1][1]),int(pt_out[0][1][0]),2], 50);

        # Verify that near by pixels are black
        self.assertLess(image_out[int(pt_out[0][1][1]-10),int(pt_out[0][1][0]-10),2],50);

        # Verify that the old pixel is now black
        self.assertLess(image_out[int(rotated_pt[0][1][1]),int(rotated_pt[0][1][0]),2],50);

    def test_draw_fluorescent_image(self):
        # Set input image
        image = cv2.cvtColor(cv2.imread(
            'image_processing/fluorescence_image_example.png'),cv2.COLOR_BGR2RGB)
        fluorescence_image_points_on_line_pix = [
                [ [ 205, 712], [ 217, 748] ], # A few points on line 0 (v)
                [ [ 325, 656], [ 328, 690] ], # A few points on line 1 (v)
                [ [ 392, 618], [ 418, 645] ], # ...
                [ [ 826, 360], [ 848, 385] ],
                [ [ 973, 280], [ 984, 310] ],
                [ [1009, 258], [1022, 294] ],
                [ [1068, 194], [1080, 233] ], # A few points on line 6 (h)
                ]
        photobleach_line_group = 'vvvhhhh'
        
        draw_fluorescent_image(
            image,
            fluorescence_image_points_on_line_pix,
            photobleach_line_group)
        
if __name__ == '__main__':
    unittest.main()
