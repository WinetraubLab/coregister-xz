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
            self.recomended_center_pix = np.array([0,0])
        else:
            raise ValueError("Invalid initialization method: %s" % method)
            
        self._check_u_v_consistency_assumptions()
    
    def _fit_from_points_on_photobleach_lines(self, 
        fluorescence_image_points_on_line_pix, photobleach_line_position_mm, photobleach_line_group):
        """ This function initialize FitPlane by points on photobleach lines.
        It sets self values of u,v,h and recomended_center_pix point in the fluorescence image c=(cu,cv)"""
        
        if (len(fluorescence_image_points_on_line_pix) != len(photobleach_line_position_mm) or
            len(fluorescence_image_points_on_line_pix) != len(photobleach_line_group)):
            raise ValueError("Number of lines should be the same between " + 
                "fluorescence_image_points_on_line_pix, photobleach_line_position_mm, photobleach_line_group")
        
        # Solve x,y first
        self._fit_from_photobleach_lines_xy(
            fluorescence_image_points_on_line_pix, 
            photobleach_line_position_mm, 
            photobleach_line_group)
            
        # Make sure u has no z component. It will help make things standard
        self._fit_from_photobleach_lines_z_from_no_shear_equal_size()
        
        # Fix z component
        self.h[2] = 0
        
        # Check
        self._check_u_v_consistency_assumptions()
        
        # Find recomended_center according to this logic:
        # c_u - according to the location that norm hits the plane
        # c_v - center of the fluorescence_image_points_on_line_pix
        v_coordinates = np.array([item[1] for sublist in fluorescence_image_points_on_line_pix for item in sublist])
        c_v = np.mean(v_coordinates)
        o = self.get_uv_from_xyz([0,0,0])
        c_u = o[0]
        self.recomended_center_pix = np.array([c_u, c_v])
        
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

    def _check_u_v_consistency_assumptions(self, skip_value_cheks=False):
        """ Check u,v assumptions """
        
        # Skip
        if skip_value_cheks:
            return
    
        # Check u and v are orthogonal and have the same norm
        if not (np.abs(self.u_norm_mm() - self.v_norm_mm())/self.v_norm_mm() < 0.05):
            raise ValueError('u and v should have the same norm')
        if not (np.dot(self.u,self.v)/(self.u_norm_mm()*self.v_norm_mm()) < 0.05):
            raise ValueError('u must be orthogonal to v')
        
        # Check that u vec is more or less in the x-y plane
        min_ratio = 0.1 # corresponding <6 degrees
        if not ( abs(self.u[2]) < np.linalg.norm(self.u[:2])*min_ratio):
            raise ValueError(
                'Make sure that tissue surface is parallel to x axis (<%.2f slope), angle is too steep right now'
                % min_ratio)

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
        a,b,c are unitless and d has units of mm """
        normal_vec = self.normal_direction()
        a, b, c = normal_vec
        d = -np.dot(normal_vec, self.h)
        
        return a,b,c,d
        
    def get_xyz_from_uv(self, point_pix):
        """ Get the 3D physical coordinates of a specific pixel in the image [u_pix, v_pix] """
        u_pix = point_pix[0]
        v_pix = point_pix[1]
        return (self.u*u_pix + self.v*v_pix + self.h)
    
    def get_uv_from_xyz(self, point_mm):
        """ Get the u,v coordinates on an image from a point in space, if point is outside the plane, return the u,v of the closest point """
        
        # Assuming u is orthogonal to v (as it shuld) for this function to work
        self._check_u_v_consistency_assumptions()
        
        u_hat = self.u_direction()
        u_norm = self.u_norm_mm()
        u_pix = np.dot(point_mm-self.h,u_hat)/u_norm
        
        v_hat = self.v_direction()
        v_norm = self.v_norm_mm()
        v_pix = np.dot(point_mm-self.h,v_hat)/v_norm
        
        return np.array([u_pix, v_pix])
        
    def get_xy_projection(self, min_x_mm=None, max_x_mm=None, min_y_mm=None, max_y_mm=None):
        """ When lookin at the pattern from above, return two points from the projected line 
        The line would go from point1 --> point2 where u value increases.
        
        USAGE:
            pt1, pt2 = get_xy_projection()
            pt1, pt2 = get_xy_projection(min_x_mm = 0, max_x_mm = 10, min_y_mm = 0, max_y_mm = 10)
            
        INPUTS:
            If none of the optional inputs are defined then line will be (u,v): (0,0) --> (c_u,c_v)
            If min_x_mm, max_x_mm are defined pt1[0] = min_x_mm, pt2[0] = max_x_mm. 
            If min_y_mm, max_y_mm are defined pt1[1] = min_y_mm, pt2[1] = max_y_mm. 
            If both sets of x and y are defined, we will use the outmost inclusive set
        """
        # Get the points on the plane that satisfy the x condition
        no_x_limit = min_x_mm is None or max_x_mm is None
        if no_x_limit:
            # No clear user limits
            pt1_u_x = np.Inf
            pt2_u_x = -np.Inf
        else:
            # We need to find where min_x_mm, max_x_mm are on the plane.
            # To do so, we get the equation ax+by+cz+d=0, and set x to the limits, and z to 0 to find y.
            a,b,c,d = self.plane_equation()
            min_x_y_mm = -(d+a*min_x_mm)/b
            max_x_y_mm = -(d+a*max_x_mm)/b
            
            # Find u,v on that plane
            tmp1 = self.get_uv_from_xyz([min_x_mm, min_x_y_mm, 0])
            tmp2 = self.get_uv_from_xyz([max_x_mm, max_x_y_mm, 0])
            pt1_u_x, pt2_u_x = min(tmp1[0],tmp2[0]), max(tmp1[0],tmp2[0])
        
        # Get the points on the plane that satisfy the y condition
        no_y_limit = min_y_mm is None or max_y_mm is None
        if no_y_limit:
            # No clear user limits
            pt1_u_y = np.Inf
            pt2_u_y = -np.Inf
        else:
            # We need to find where min_y_mm, max_y_mm are on the plane.
            # To do so, we get the equation ax+by+cz+d=0, and set y to the limits, and z to 0 to find x.
            a,b,c,d = self.plane_equation()
            min_y_x_mm = -(d+b*min_y_mm)/a
            max_y_x_mm = -(d+b*max_y_mm)/a
            
            # Find u,v on that plane
            tmp1 = self.get_uv_from_xyz([min_y_x_mm, min_y_mm, 0])
            tmp2 = self.get_uv_from_xyz([max_y_x_mm, max_y_mm, 0])
            pt1_u_y, pt2_u_y = min(tmp1[0],tmp2[0]), max(tmp1[0],tmp2[0]) 

        if no_x_limit and no_y_limit:
            # No limits found, use default values 
            pt1_u, pt2_u = min(0, self.recomended_center_pix[0]), max(0, self.recomended_center_pix[0])
        else:        
            # Aggregate all points to find the maximum bounds
            pt1_u = min(pt1_u_x,pt1_u_y)
            pt2_u = max(pt2_u_x,pt2_u_y)
        pt12_v = self.recomended_center_pix[1]
        
        # Figure out u,v on the plane that the points correspond to
        pt1 = self.get_xyz_from_uv([pt1_u, pt12_v])
        pt2 = self.get_xyz_from_uv([pt2_u, pt12_v])
        
        return (pt1[:2],pt2[:2])
