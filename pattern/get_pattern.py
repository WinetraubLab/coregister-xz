import numpy as np

""" Pattern parameters """
# Line Positions
v_line_bias_mm = 500e-3
h_line_bias_mm = -650e-3
base_mm = 90e-3
# Unitless, units will be added later in the code. Make sore it is sorted.
v_lines_spacing = [-1, 0, 1]
# Unitless, units will be added later in the code. Make sure it is sorted.
h_lines_spacing = [-3, 0, 1, 3]


def get_number_of_lines():
    return (len(v_lines_spacing) + len(h_lines_spacing)))

def get_vh_photobleach_line_position_mm():
    v_photobleach_line_position_mm = np.array(v_lines_spacing) * base_mm + v_line_bias_mm
    h_photobleach_line_position_mm = np.array(h_lines_spacing) * base_mm + h_line_bias_mm

    # Check sorting otherwise the relationship between cut position and the graphs in the last cell is lost.
    if not np.array_equal(v_photobleach_line_position_mm, np.sort(v_photobleach_line_position_mm)):
      raise ValueError('v_lines_spacing should be sorted')
    if not np.array_equal(h_photobleach_line_position_mm, np.sort(h_photobleach_line_position_mm)):
      raise ValueError('h_lines_spacing should be sorted')

    return v_photobleach_line_position_mm, h_photobleach_line_position_mm


""" This function returns the pattern as seen by fluorescent image (left to right)
    INPUTS:
        cut_position - can be 1,2,3,4
    OUTPUTS:
        photobleach_line_position_mm
        photobleach_line_group
"""
def get_pattern_on_fluorescent_image(cut_position=2):

    v_photobleach_line_position_mm, h_photobleach_line_position_mm = get_vh_photobleach_line_position_mm()

    # Combine v lines and h lines according to order
    def combine(v_photobleach_line_position_mm, h_photobleach_line_position_mm,
                order='vh',flip_v=False, flip_h=False):
      if flip_v:
        v = np.flipud(v_photobleach_line_position_mm)
      else:
        v = v_photobleach_line_position_mm

      if flip_h:
        h = np.flipud(h_photobleach_line_position_mm)
      else:
        h = h_photobleach_line_position_mm

      if order == 'hv':
        photobleach_line_position_mm = np.concatenate((h, v), axis=0)
        photobleach_line_group = (['h'] * len(h)) + (['v'] * len(v))
      else:
        photobleach_line_position_mm = np.concatenate((v, h), axis=0)
        photobleach_line_group = (['v'] * len(v)) + (['h'] * len(h))

      return photobleach_line_position_mm, photobleach_line_group

    if cut_position == 1:
      photobleach_line_position_mm, photobleach_line_group = combine(
          v_photobleach_line_position_mm, h_photobleach_line_position_mm,
          order='vh', flip_v=False, flip_h=True)
    elif cut_position == 2:
      photobleach_line_position_mm, photobleach_line_group = combine(
          v_photobleach_line_position_mm, h_photobleach_line_position_mm,
          order='hv', flip_h=False, flip_v=False)
    elif cut_position == 3:
      photobleach_line_position_mm, photobleach_line_group = combine(
          v_photobleach_line_position_mm, h_photobleach_line_position_mm,
          order='hv', flip_h=True, flip_v=False)
    elif cut_position == 4:
      photobleach_line_position_mm, photobleach_line_group = combine(
          v_photobleach_line_position_mm, h_photobleach_line_position_mm,
          order='vh', flip_v=False, flip_h=False)

    return photobleach_line_position_mm, photobleach_line_group
