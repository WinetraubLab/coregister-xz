import numpy as np
from scipy.optimize import minimize
import cv2
import matplotlib.pyplot as plt


def draw_fluorescent_image(image,fluorescence_image_points_on_line_pix,photobleach_line_group):
  # Plot image
  plt.imshow(image)

  # Plot lines
  for points, group in zip(fluorescence_image_points_on_line_pix, photobleach_line_group):
    if group == 'v' or group == 'V':
      color = 'red'
    else:
      color = 'blue'

    points = np.array(points)
    plt.plot(points[:,0], points[:,1], color=color, linewidth=1)

  plt.axis('off')  # Hide the axis
  plt.show()
