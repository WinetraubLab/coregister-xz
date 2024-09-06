import numpy as np
import matplotlib.pyplot as plt
from plane.fit_plane import FitPlane

def plot_fit_plane_xy(
    fp, #FitPlane or array of FitPlanes
    v_lines_mm, # Position of vertical lines
    h_lines_mm, # Position of horizontal lines
    oct_scan_size_mm=0.5, # Size of the OCT scan around the center
    plot_bound_mm=2.5, # How big to plot
    reverse_plot=False, # Set to true if the flourecence image is reversed
    ax=None, # Set to axs if plotting in a subplot is needed, use None otherwise
    ):
    """ Plot the fit plane from above (xy projection) """
    
    # Input check
    if not isinstance(fp, list):
        fp = [fp]

    if ax is None:
        _, ax = plt.subplots()
        show_at_end = True
    else:
        show_at_end = False
    
    # Plot photobleach lines pattern
    for v_line in v_lines_mm:
        ax.axvline(x=v_line, color='r', linestyle='-')
    for h_line in h_lines_mm:
        ax.axhline(y=h_line, color='b', linestyle='-')
    
    # Plot OCT Scan
    square_x = [-oct_scan_size_mm/2, oct_scan_size_mm/2, oct_scan_size_mm/2, -oct_scan_size_mm/2, -oct_scan_size_mm/2]
    square_y = [-oct_scan_size_mm/2, -oct_scan_size_mm/2, oct_scan_size_mm/2, oct_scan_size_mm/2, -oct_scan_size_mm/2]
    ax.plot(square_x, square_y, color='k', linestyle=':')
    
    # Draw the cut planes (arrows)
    for fp_instance in fp:
        pt1,pt2 = fp_instance.get_fit_plane_xy_projection(
            min_x_mm = min(v_lines_mm)-0.1,
            max_x_mm = max(v_lines_mm)+0.1,
            min_y_mm = min(h_lines_mm)-0.1,
            max_y_mm = max(h_lines_mm)+0.1,
            )
        d = pt2-pt1
        ax.arrow(pt1[0], pt1[1], d[0], d[1], color='k', head_width=0.1, head_length=0.1)
    
    # Set titles, axis etc
    ax.set_xlabel('X[mm]')
    ax.set_ylabel('Y[mm]')
    ax.grid(True)
    ax.axis('square')
    ax.set_xlim(-plot_bound_mm, plot_bound_mm)
    ax.invert_yaxis() # Plot using ij notation instead of xy
    if reverse_plot:
        ax.invert_yaxis()
        ax.invert_xaxis()

    if show_at_end:
        plt.show()
