import numpy as np
from scipy.spatial import cKDTree  # type: ignore
from matplotlib import path

def create_nn_interpFunc(x_coords, y_coords, grid_size, k_nn=1, 
                         polygon_csv=None, extents=None,
                        ):
    """
    Creates a nearest-neighbor interpolator function for a fixed grid.
    
    Parameters:
    - coords: List of (x, y) coordinate pairs (length L).
    - grid_size: Desired grid cell size.
    - k_nn: Number of nearest neighbors to use (default is 1).
    - polygon_csv: Path to csv file containing polygon vertices (columns: x, y). 
                   If provided, values outside the polygon will be masked as NaN.
    - extents: Manual x- and y-extent limits for the interpolation grid.
                   Example: extents=(10, 50, 20, 80).
                   If None, this defaults to the min/max of the input coords.
    - tree_inds, weights, inside_mask: Optional precomputed values to speed up repeated calls.

    Returns:
    - A function that can be used to interpolate new z values efficiently.
    """

    # Convert input coordinates to NumPy array
    x_coords, y_coords = np.asarray(x_coords), np.asarray(y_coords)

    # Set extents to full mesh extent if smaller window not provided
    if extents is None:
        extents = (x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max())

    ext_inds = np.where((extents[0] <= x_coords) & (x_coords <= extents[1]) &
                        (extents[2] <= y_coords) & (y_coords <= extents[3]))[0]
    
    # Filter coordinates to those within extents
    x_coords, y_coords = x_coords[ext_inds], y_coords[ext_inds]

    # Derive grid limits from coords, which have already been filtered to the desired extent
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()

    # Create the grid
    x_grid = np.arange(x_min, x_max + grid_size, grid_size)
    y_grid = np.arange(y_min, y_max + grid_size, grid_size)
    X, Y = np.meshgrid(x_grid, y_grid, indexing="xy")
    grid_points = np.column_stack([X.ravel(), Y.ravel()])  # Flatten grid

    #if tree_inds is None and inside_mask is None:
    # Build KDTree and find nearest-neighbor indices
    tree = cKDTree([(x, y) for x, y in zip(x_coords, y_coords)])
    distances, tree_inds = tree.query(grid_points, k=k_nn, workers=-1)  # Fast parallel query

    # ---- Polygon mask (optional) ---- #
    if polygon_csv is None:
        inside_mask = np.ones(X.shape, dtype=bool)
    else:
        # Load polygon vertices (supports headers)
        poly_data = np.loadtxt(polygon_csv, delimiter=",")
        poly_x, poly_y = poly_data[:, 0], poly_data[:, 1]

        poly_path = path.Path(np.column_stack((poly_x, poly_y)))
        inside_mask = poly_path.contains_points(grid_points).reshape(X.shape)
        

    def nn_interpolator(z_values):
        """
        Interpolates z-values using k-nearest neighbors.

        Parameters:
        - z_values: Array of shape (L,) corresponding to coords.

        Returns:
        - Z: 2D array of interpolated values matching the grid shape.
        """

        # Filter z_values to match extents
        z_values = np.asarray(z_values)[ext_inds]
        assert len(z_values) == len(x_coords), "Mismatch in z-values length"

        if k_nn == 1:
            Z_interp = z_values[tree_inds]
        else:
            # Weight by inverse distance, avoiding divide-by-zero
            with np.errstate(divide='ignore'):
                weights = 1./distances
            weights[np.isinf(weights)] = 1e10  # handle exact matches
            weights /= np.sum(weights, axis=1, keepdims=True)
            Z_interp = np.sum(z_values[tree_inds]*weights, axis=1)

        Z = Z_interp.reshape(Y.shape)

        # Apply mask
        Z[~inside_mask] = np.nan

        # Flip vertically to match typical image orientation
        Z = np.flipud(Z)

        return Z

    return nn_interpolator
