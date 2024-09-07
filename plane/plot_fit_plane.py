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

def plot_fit_plane_uv(
    fp, #FitPlane or array of FitPlanes
    v_lines_mm, # Position of vertical lines
    h_lines_mm, # Position of horizontal lines
    v_range, # tuple (min_v, max_v) in pixels
    ax=None, # Set to axs if plotting in a subplot is needed, use None otherwise
    ):
    """ Plot the fit plane from above (uv projection) """

    # Input checks
    if ax is None:
        _, ax = plt.subplots()
        show_at_end = True
    else:
        show_at_end = False

    # Resize
    fig = ax.figure
    fig.set_size_inches(10, 2)

    # Intercept satisfies a*u+b*v+c=0
    def plot_intercept(a, b, c, color, label):
        # Compute points u,v
        v1, v2 = v_range
        u1 = -(c+b*v1)/a
        u2 = -(c+b*v2)/a

        # Plot
        ax.plot([u1, u2],[v1, v2],color=color, linestyle='-', label=label)

    # Loop over photobleach lines v and h
    for v_line in v_lines_mm:
        a,b,c = fp.get_v_line_fit_plane_intercept(v_line)
        plot_intercept(a,b,c,'r','')
    plot_intercept(a,b,c,'r','v')
    for h_line in h_lines_mm:
        a,b,c = fp.get_h_line_fit_plane_intercept(h_line)
        plot_intercept(a,b,c,'b','')
    plot_intercept(a,b,c,'b','h')
  
    # Set titles, axis etc
    ax.set_xlabel('U[pix]')
    ax.set_ylabel('V[pix]')
    ax.grid(True)
    ax.set_aspect('equal')
    ax.legend(loc='upper center', bbox_to_anchor=(0.8, -0.2), ncol=2)

    if show_at_end:
        plt.show()

