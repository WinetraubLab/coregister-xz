import numpy as np

""" This class defines a photobleach pattern and provide basic functionality """
class Pattern:
    def __init__ (
        self,
        v_line_bias_mm = 500e-3,
        h_line_bias_mm = -650e-3,
        base_mm = 90e-3,
        v_line_positions = [-1, 0, 1], # Unitless, units will be added later in the code. Make sore it is sorted.
        h_line_positions = [-3, 0, 1, 3],
        ):

        # Input check
        # Check sorting otherwise the relationship between cut position and the graphs in the last cell is lost.
        if not np.array_equal(v_line_positions, np.sort(v_line_positions)):
          raise ValueError('v_line_positions should be sorted')
        if not np.array_equal(h_line_positions, np.sort(h_line_positions)):
          raise ValueError('h_line_positions should be sorted')
        
        # Assign
        self.v_line_bias_mm = v_line_bias_mm
        self.h_line_bias_mm = h_line_bias_mm
        self.base_mm = base_mm
        self.v_line_positions = v_line_positions
        self.h_line_positions = h_line_positions

        # Compute line positions in mm
        self.v_line_positions_mm = np.array(v_line_positions) * base_mm + v_line_bias_mm
        self.h_line_positions_mm = np.array(h_line_positions) * base_mm + h_line_bias_mm

    def get_number_of_lines(self):
        return (len(self.v_line_positions) + len(self.h_line_positions))


    """ This function returns the pattern as seen by fluorescent image (left to right)
        INPUTS:
            cut_position - can be 1,2,3,4
        OUTPUTS:
            photobleach_line_position_mm
            photobleach_line_group
    """
    def project_pattern_onto_plane(self,cut_position=2):

        # Combine v lines and h lines according to order
        def combine(order='vh',flip_v=False, flip_h=False):
          if flip_v:
            v = np.flipud(self.v_line_positions_mm)
          else:
            v = self.v_line_positions_mm

          if flip_h:
            h = np.flipud(self.h_line_positions_mm)
          else:
            h = self.h_line_positions_mm

          if order == 'hv':
            photobleach_line_position_mm = np.concatenate((h, v), axis=0)
            photobleach_line_group = (['h'] * len(h)) + (['v'] * len(v))
          else:
            photobleach_line_position_mm = np.concatenate((v, h), axis=0)
            photobleach_line_group = (['v'] * len(v)) + (['h'] * len(h))

          return photobleach_line_position_mm, photobleach_line_group

        if cut_position == 1:
          photobleach_line_position_mm, photobleach_line_group = combine(
              order='vh', flip_v=False, flip_h=True)
        elif cut_position == 2:
          photobleach_line_position_mm, photobleach_line_group = combine(
              order='hv', flip_h=False, flip_v=False)
        elif cut_position == 3:
          photobleach_line_position_mm, photobleach_line_group = combine(
              order='hv', flip_h=True, flip_v=False)
        elif cut_position == 4:
          photobleach_line_position_mm, photobleach_line_group = combine(
              order='vh', flip_v=False, flip_h=False)

        return photobleach_line_position_mm, photobleach_line_group
