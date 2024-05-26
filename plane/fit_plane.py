import numpy as np

class FitPlane:
    
    def __init__(self, *args, method="points on photobleach lines") :
        """
        Constructor for the FitPlane class based on method

        If method="points on photobleach lines", Inputs:
            fluorescence_image_points_on_line_pix: For each of the photobleach lines, 
                find at least two points in the fluorescence image. Mark the coordinates as pixel values
                l1 = [[x1,y1],[x2,y2],...] and create an array of those [l1,l2,...,ln]
            photobleach_line_position_mm: an array defining the position (in mm) of each of the photobleach line positions 
            photobleach_line_group: an array defining each line is a horizontal or vertical line ['h','v',...] etc
        If method="u,v,h directly", Inputs:
            u,v,h in mm
        """
        if method == "points on photobleach lines":
            self._fit_from_points_on_photobleach_lines(args[0],args[1],args[2])
        elif method == "u,v,h directly":
            self.u = np.array(args[0])
            self.v = np.array(args[1])
            self.h = np.array(args[2])
        else:
            raise ValueError("Invalid initialization method: %s" % method)
    
    def _fit_from_points_on_photobleach_lines(self, 
        fluorescence_image_points_on_line_pix, photobleach_line_position_mm, photobleach_line_group,
        override_value_cheks=False):
        
        # Solve x,y first
        self._fit_from_photobleach_lines_xy(
            fluorescence_image_points_on_line_pix, 
            photobleach_line_position_mm, 
            photobleach_line_group)
            
        # Make sure u has no z component. It will help make things standard
        self._fit_from_photobleach_lines_z_from_no_shear_equal_size()
        
        # Fix z component
        self.h[2] = 0
        
        # Check that u vec is more or less in the x-y plane
        min_ratio = 0.1 # corresponding <6 degrees
        if (not ( abs(self.u[2]) < np.linalg.norm(self.u[:2])*min_ratio) and 
            not (override_value_cheks)):
            raise ValueError(
                'Make sure that tissue surface is parallel to x axis (<%.2f slope), angle is too steep right now'
                % min_ratio)
        
    def _fit_from_photobleach_lines_xy(self, 
        fluorescence_image_points_on_line_pix, photobleach_line_position_mm, photobleach_line_group
        ):
        """ First part estimates x-y part of u,v,h using least squares matrix """
        
        # Generate least square matrix
        def gen_row (
            ln_pts, # Points on the specific line
            ln_id_group, # Can be 'h' or 'v'
            ln_id_pos,   # The physical position in mm
            ):
            y_row = []
            A_row = []
            for ln_pt in ln_pts:
                y.append(ln_id_pos) # least square y is the line position [mm]
                if ln_id_group == 'v':
                    A_row.append([ln_pt[0], 0, ln_pt[1], 0, 1, 0])
                else:
                    A_row.append([0, ln_pt[0], 0, ln_pt[1], 0, 1])

            A_row = np.array(A_row)
            y_row = np.array(y_row)

            return (A_row,y_row)
        A = []
        y = []
        for index, _ in enumerate(photobleach_line_group):
            A_row, y_row = gen_row(fluorescence_image_points_on_line_pix[index],
                photobleach_line_group[index],photobleach_line_position_mm[index])
            A.append(A_row)
            y.append(y_row)
        A = np.vstack(A)
        y = np.hstack(y)
        
        # Solve the least square problem
        x, residuals, rank, s = np.linalg.lstsq(A, y, rcond=None)

        # Output vectors
        self.u = np.array([x[0], x[1], np.nan])
        self.v = np.array([x[2], x[3], np.nan])
        self.h = np.array([x[4], x[5], np.nan])
        
    def _fit_from_photobleach_lines_z_from_no_shear_equal_size(self):
        # Estimate z by solving equation A and equation B
        u_x = self.u[0]
        u_y = self.u[1]
        v_x = self.v[0]
        v_y = self.v[1]
        def eq_A(u_x,u_y,v_x,v_y):
            return u_x*v_x+u_y*v_y
        def eq_B(u_x,u_y,v_x,v_y):
            return u_x**2-v_x**2+u_y**2-v_y**2
        def eq_vz(u_x,u_y,v_x,v_y):
            A = eq_A(u_x,u_y,v_x,v_y)
            B = eq_B(u_x,u_y,v_x,v_y)
            return 1/np.sqrt(2)*np.sqrt(B+np.sqrt(B**2-4*A**2))
            
        # Estimate z component
        v_z = eq_vz(u_x,u_y,v_x,v_y)
        self.u[2] = -eq_A(u_x,u_y,v_x,v_y)/v_z
        self.v[2] = v_z
        
        # Check consistency assumptions
        assert(np.abs(self.u_norm_mm() - self.v_norm_mm())/self.v_norm_mm() < 0.05)
        assert(np.dot(self.u,self.v)/(self.u_norm_mm()*self.v_norm_mm()) < 0.05)

    def u_norm_mm(self):
        """ Return the size of pixel u in mm """
        return np.linalg.norm(self.u)
    
    def v_norm_mm(self):
        """ Return the size of pixel v in mm """
        return np.linalg.norm(self.v)
    
    def u_direction(self):
        """ Return a unit vector in the direction of u """
        return self.u / self.u_norm_mm()
        
    def v_direction(self):
        """ Return a unit vector in the direction of v """
        return self.v / self.v_norm_mm()
        
    def normal_direction(self):
        """ Return a unit vector in the direction of the normal """
        return np.cross(self.u_direction(), self.v_direction())
        
    def plane_equation(self):
        """ Convert u,v,h to a plane equation ax+by+cz-d=0.
        a,b,c are unitless and d has units of mm"""
        normal_vec = self.normal_direction()
        a, b, c = normal_vec
        d = -np.dot(normal_vec, self.h)
        
        return a,b,c,d
        
    
