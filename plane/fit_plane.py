import numpy as np

class FitPlane:
    
    def __init__(self, 
        fluorescence_image_points_on_line_pix, photobleach_line_position_mm,photobleach_line_group
        ):
        """
        Constructor for the FitPlane class.

        :param fluorescence_image_points_on_line_pix: For each of the photobleach lines, 
            find at least two points in the fluorescence image. Mark the coordinates as pixel values
            l1 = [[x1,y1],[x2,y2],...] and create an array of those [l1,l2,...,ln]
        :param photobleach_line_position_mm: an array defining the position (in mm) of each of the photobleach line positions 
        :param photobleach_line_group: an array defining each line is a horizontal or vertical line ['h','v',...] etc
        """
        
        self._fit_from_photobleach_lines(
            fluorescence_image_points_on_line_pix, photobleach_line_position_mm,photobleach_line_group)
    
    def _fit_from_photobleach_lines(self, 
        fluorescence_image_points_on_line_pix, photobleach_line_position_mm,photobleach_line_group
        ):
        
        self.u = np.array([0.707,-0.707,0])
        self.v = np.array([0,0,1])
        self.h = np.array([0,0,0])

    def u_norm(self):
        """ Return the size of pixel u in mm """
        return np.linalg.norm(self.u)
    
    def v_norm(self):
        """ Return the size of pixel v in mm """
        return np.linalg.norm(self.v)
    
    def u_direction(self):
        """ Return a unit vector in the direction of u """
        return self.u / self.u_norm()
        
    def v_direction(self):
        """ Return a unit vector in the direction of v """
        return self.v / self.v_norm()
