import numpy as np
import matplotlib.pyplot as plt
from plane.fit_plane import FitPlane

def plot_fit_plane(
    fp, #FitPlane or array of FitPlanes
    vLines_mm, # Position of vertical lines
    hLines_mm, # Position of horizontal lines
    oct_scan_size_mm=0.5, # Size of the OCT scan around the center
    plot_bound_mm=2.5, # How big to plot
    reverse_plot=False, # Set to true if the flourecence image is reversed
    ):
    """ Plot the fit plane from above (xy projection) """
    
    # Input check
    if not isinstance(fp, list):
        fp = [fp]
    
    # Plot photobleach lines pattern
    for vline in vLines_mm:
        plt.axvline(x=vline, color='r', linestyle='-')
    for hline in hLines_mm:
        plt.axhline(y=hline, color='b', linestyle='-')
    
    # Plot OCT Scan
    square_x = [-oct_scan_size_mm/2, oct_scan_size_mm/2, oct_scan_size_mm/2, -oct_scan_size_mm/2, -oct_scan_size_mm/2]
    square_y = [-oct_scan_size_mm/2, -oct_scan_size_mm/2, oct_scan_size_mm/2, oct_scan_size_mm/2, -oct_scan_size_mm/2]
    plt.plot(square_x, square_y, color='k', linestyle=':')
    
    # Draw the planes
    for fp_instance in fp:
        pt1,pt2 = fp_instance.get_xy_projection(
            min_x_mm = min(vLines_mm)-0.1,
            max_x_mm = max(vLines_mm)+0.1,
            min_y_mm = min(hLines_mm)-0.1,
            max_y_mm = max(hLines_mm)+0.1,
            )
        d = pt2-pt1
        plt.arrow(pt1[0], pt1[1], d[0], d[1], color='k', head_width=0.1, head_length=0.1)
    
    # Set titles, axis etc
    plt.xlabel('X[mm]')
    plt.ylabel('Y[mm]')
    plt.grid(True)
    plt.axis('equal')
    plt.xlim(-plot_bound_mm, plot_bound_mm)
    #plt.ylim(-plot_bound_mm, plot_bound_mm)
    if reverse_plot:
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()

    plt.show()
