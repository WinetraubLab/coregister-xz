import numpy as np
from scipy.optimize import minimize
import cv2
from plane.fit_plane import FitPlane


""" This file contains utility functions to pre-process fluorescent images
with barcodes to make sure they meet _check_u_v_consistency_assumptions set by fitplane """


""" Rotate image to meet requirements
    For explanation about this function see: https://docs.google.com/document/d/1-JUh07n5B2JEqUmg6tkGTbQ7zrT5hNYcfUa5oMXOGm4/edit?tab=t.0#bookmark=id.ab2crgzg4856
    INPUTS:
    :fluorescence_image_points_on_line_pix: List of points. Dimensions are
        [line number][point number on the line][0=u,1=v]
    :photobleach_line_group: an array defining each line is a horizontal or vertical line ['h','v',...]
    :fluorescence_image: cv object with the image

    OUTPUT: Same as input but rotated to maximize u-v consistency

    """
def rotate_image_to_meet_consistency_assumptions(
    fluorescence_image_points_on_line_pix,
    photobleach_line_group,
    fluorescence_image):

    # Get the angle
    theta_rad = _rotate_points_to_optimize_uv(
        fluorescence_image_points_on_line_pix,
        photobleach_line_group)

    # Rotate points along the center of the image
    fluorescence_image_points_on_line_pix_out = _rotate_pt_list(
        fluorescence_image_points_on_line_pix, theta_rad,
        center_x = (fluorescence_image.shape[1]-1)/2,
        center_y = (fluorescence_image.shape[0]-1)/2,
        )

    # Rotate image - calculate the rotation matrix
    M = cv2.getRotationMatrix2D(((fluorescence_image.shape[1]-1)/2, (fluorescence_image.shape[0]-1)/2),
                                np.degrees(theta_rad), 1)

    # Apply the rotation
    fluorescence_image_out = cv2.warpAffine(
        fluorescence_image, M, (fluorescence_image.shape[1], fluorescence_image.shape[0]),
        flags=cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT)

    return fluorescence_image_points_on_line_pix_out, fluorescence_image_out


# Optimize rotation angle, output is the optimum
def _rotate_points_to_optimize_uv(
    fluorescence_image_points_on_line_pix,
    photobleach_line_group):

    def error_fun(theta_rad, is_print=False):
        theta_rad = theta_rad[0]
            
        pt_rot = _rotate_pt_list(
            fluorescence_image_points_on_line_pix, theta_rad)

        fp = FitPlane.from_fitting_points_on_photobleach_lines(
            pt_rot, range(len(pt_rot)),
            photobleach_line_group,
            skip_uv_value_cheks=True,
            )

        # Value to minimize
        e = (fp.v[0]**2+fp.v[1]**2)/(fp.u[0]**2+fp.u[1]**2)

        if is_print:
            print(f"Theta: {theta_rad}, e={e}")
            print(fp.u_direction())
            print(fp.v_direction())

        return e

    # Optimize
    initial_guess = 0
    result = minimize(
        error_fun, initial_guess,
        bounds=[(-np.pi, np.pi)])
    theta_out = result.x[0]

    # Extract the result
    return theta_out
        
# Rotate points acording to angle
def _rotate_pt_list(pt_list, theta_rad, center_x=0, center_y=0):
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    m = np.array([[c,s],[-s,c]])

    p0 = np.array([center_x, center_y])
        
    new_pt_list=[]
    for line_points in pt_list:
        new_line_points = []
        for point in line_points:
            # Rotate vector using p0 as the center
            v = m.dot(np.array(point)-p0)+p0

            # Store result
            new_line_points.append(v.tolist())
        new_pt_list.append(new_line_points)

    return new_pt_list
