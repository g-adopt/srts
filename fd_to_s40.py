""" Adopted from the script Tom New has provided
"""
# %%
import numpy as np
from gdrift.profile import SplineProfile
import pyvista as pv
from pathlib import Path
import gdrift
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmcrameri.cm as cmc
from matplotlib.colors import TwoSlopeNorm, Normalize


def __main__(should_plot_layer: bool = False):
    # Dimensional constants for conversion and model
    dimensional_constants = get_dimensional_constants()

    # Path to the pvtu file -- this is a gadopt model
    pvtu_file = Path(
        "output_0.pvtu")

    # path to models directory in S40RTS code
    name = "0Ma_C52_8e8_Kappa"
    output_dir = Path(f"./geodyn/{name}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # load the model
    model = load_pvtu_points(pvtu_file)
    model, avg_temperature_spline = convert_model(model)

    # From this point on we are interpolating into a depth profile consitent with
    # srts example model
    # create array of interpolation depths
    # Just a rough number to make it equivalent to usual asumption of lon lat grids
    equidistanced_xyz_sphere = gdrift.fibonacci_sphere(n=180 * 360)
    _, lons, lats = gdrift.cartesian_to_geodetic(
        equidistanced_xyz_sphere[:, 0],
        equidistanced_xyz_sphere[:, 1],
        equidistanced_xyz_sphere[:, 2]
    )
    lons = np.linspace(-180, 180, dimensional_constants["n_lons"])
    lats = np.linspace(-90, 90, dimensional_constants["n_lats"])

    depth_profile = np.linspace(
        0, gdrift.R_earth - gdrift.R_cmb, dimensional_constants["nlayers"])

    depth_grid, lon_grid, lat_grid = np.meshgrid(
        depth_profile, lons, lats, indexing="ij")

    model_lat, model_lon, model_depth = gdrift.cartesian_to_geodetic(
        *[model.points[:, i] * gdrift.R_earth / dimensional_constants["r_max"] for i in range(3)])

    dists, inds = cKDTree(
        np.column_stack(
            (
                np.asarray(model_depth)/1e3,
                np.asarray(model_lon),
                np.asarray(model_lat)
            )
        )
    ).query(np.column_stack(((depth_grid).flatten()/1e3, lon_grid.flatten(), lat_grid.flatten())), k=10)

    # Compute Gaussian weights based on distances
    sigma = 1.0  # Gaussian width parameter - adjust as needed
    weights = np.exp(-0.5 * (dists / sigma)**2)
    weights /= np.sum(weights, axis=1, keepdims=True)  # normalize weights

    models_in_geodetic_coordinates = {}

    for key in ["Vs", "Vp", "Vs_elastic", "Vp_elastic"]:
        # Compute weighted average of dVs using einsum
        models_in_geodetic_coordinates[key] = np.einsum(
            'ij,ij->i',
            weights,
            model.point_data[key][inds]
        ).reshape(depth_grid.shape)

    # # save depths
    np.savetxt(output_dir / "depth_layers.dat",
               depth_profile/1e3, fmt='%04d', newline='\n')

    mapping_names = {
        "Vs": "dvs",
        "Vp": "dvp",
        "Vs_elastic": "dvs",
        "Vp_elastic": "dvp"
    }

    for key in ["Vs", "Vp", "Vs_elastic", "Vp_elastic"]:
        for ir, dpth in enumerate(depth_profile):
            # Computing the mean of each layer:
            layer_mean = np.mean(models_in_geodetic_coordinates[key][ir, :, :])
            # Converting to dlnVs/p (%):
            model_to_be_output = 100 * \
                (models_in_geodetic_coordinates[key]
                 [ir, :, :] - layer_mean) / layer_mean

            # force top and bottom layers to zero:
            if ir == 0 or ir == len(depth_profile) - 1:
                model_to_be_output = np.zeros_like(model_to_be_output)

            if should_plot_layer:
                plot_layer(lon_grid[ir, :, :], lat_grid[ir, :, :], model_to_be_output,
                           save_path=output_dir / f"{name}_{key}.{mapping_names[key]}.layer.{ir:03d}.png")

            output_file_dVs = output_dir / \
                f"{name}_{key}.{mapping_names[key]}.layer.{ir:03d}.dat"

            # output_file_dVp = output_dir / f"{name}.dvp.layer.{i+1:03d}.dat"
            with open(output_file_dVs, mode="w") as f:
                f.write(
                    "\n".join([f"{lon:+.2f} {lat:+.2f} {dat:.5f}"
                               for lon, lat, dat in zip(
                        lon_grid[ir, :, :].flatten(),
                        lat_grid[ir, :, :].flatten(),
                        model_to_be_output.flatten(),
                    )
                    ])
                )


def get_dimensional_constants():
    return {
        "T_0": 3700,  # Used for dimensionalising temperature T_nd * T_0 + T_1
        "T_1": 300,  # See above
        "r_max": 2.208,  # Radius of the Earth in nd
        "r_min": 1.208,  # Radius of the CMB in nd
        "nlayers": 65,
        "n_lons": 721,
        "n_lats": 361,
        "depth_res": 50e3,  # in meters
    }


def load_pvtu_points(pvtu_file):
    """ Loads the pvtu file and returns the model """
    # Load the pvtu file
    model = pv.read(pvtu_file).clean()

    needed_arrays = ["FullTemperature_CG", "Temperature_Deviation_CG"]

    # drop unneeded arrays
    for array_name in model.point_data.keys():
        if array_name not in needed_arrays:
            model.point_data.pop(array_name)

    dimensional_constants = get_dimensional_constants()

    # calculate T and T_av, dropping arrays after they become unneeded
    model.point_data["T"] = model["FullTemperature_CG"] * \
        dimensional_constants["T_0"] + dimensional_constants["T_1"]
    model.point_data["dT"] = model["Temperature_Deviation_CG"] * \
        dimensional_constants["T_0"]
    model.point_data["T_av"] = model["T"] - model["dT"]

    # depth in meters
    model.point_data["depth"] = gdrift.R_earth - np.linalg.norm(
        model.points * gdrift.R_earth / dimensional_constants["r_max"], axis=1)

    return model


def convert_model(model):
    dimensional_constants = get_dimensional_constants()
    # This is the depth profile for which we generate the thermodynamic model
    depth_profile = np.linspace(
        0, gdrift.R_earth - gdrift.R_cmb, dimensional_constants["nlayers"])
    # initialise thermodynamic model
    slb_pyrolite = gdrift.ThermodynamicModel(
        "SLB_16", "pyrolite", depths=depth_profile,
    )

    # Note: This is a very weird way of computing the average temperature profile,
    # along which we regularise the thermodynamic table. Best would
    # be to write out the layer avarage from G-ADOPT and use that as the
    # average temperature profile.
    # Here we just simply assume that `T_av` field at any point, is the avereage
    # temperature. Which is a very valid assumption!
    dists, inds = cKDTree(model.point_data["depth"][:, np.newaxis]).query(
        depth_profile[:, np.newaxis], k=10)

    # A temperautre profile representing the mantle average temperature
    weights = 1.0 / (dists + 1e-12)  # avoid division by zero
    weights /= np.sum(weights, axis=1, keepdims=True)  # normalize weights
    avg_temperature_profile = np.einsum(
        'ij,ij->i', weights, model.point_data["T_av"][inds])

    # Average temperature profile (used in the regularisation)
    avg_temperature_spline = gdrift.SplineProfile(
        depth=depth_profile,
        value=avg_temperature_profile,
    )
    # Regularisation of the thermodynamic table
    linear_slb_pyrolite = gdrift.mineralogy.regularise_thermodynamic_table(
        slb_pyrolite,
        avg_temperature_spline,
        regular_range={"v_s": [-1.0, 0.0],
                       "v_p": [-1.0, 0.0], "rho": [-1.0, 0.0]},
    )

    # cammarano_q_model = "Q6"  # choose model from cammarano et al., 2003
    solidus_ghelichkhan = build_solidus()
    cammarano_q_model = "Q6"  # choose model from cammarano et al., 2003
    anelasticity = build_anelasticity_model(
        solidus_ghelichkhan, q_profile=cammarano_q_model)
    linear_anelastic_slb_pyrolite = gdrift.apply_anelastic_correction(
        linear_slb_pyrolite, anelasticity)

    # compute seismic velocities in P and S
    model.point_data["Vs"] = linear_anelastic_slb_pyrolite.temperature_to_vs(
        temperature=np.asarray(model["T"]), depth=np.asarray(model["depth"]))
    model.point_data["Vp"] = linear_anelastic_slb_pyrolite.temperature_to_vp(
        temperature=np.asarray(model["T"]), depth=np.asarray(model["depth"]))
    model.point_data["Vs_elastic"] = linear_slb_pyrolite.temperature_to_vs(
        temperature=np.asarray(model["T"]), depth=np.asarray(model["depth"]))
    model.point_data["Vp_elastic"] = linear_slb_pyrolite.temperature_to_vp(
        temperature=np.asarray(model["T"]), depth=np.asarray(model["depth"]))

    return model, avg_temperature_spline


def build_anelasticity_model(solidus, q_profile: str = "Q1"):

    cammarano_parameters = {
        "Q1": {
            "B": [0.5, 10],
            "g": [20, 10]
        },
        "Q2": {
            "B": [0.8, 15],
            "g": [20, 10]
        },
        "Q3": {
            "B": [1.1, 20],
            "g": [20, 10]
        },
        "Q4": {
            "B": [0.035, 2.25],
            "g": [30, 15]
        },
        "Q5": {
            "B": [0.056, 3.6],
            "g": [30, 15]
        },
        "Q6": {
            "B": [0.077, 4.95],
            "g": [30, 15]
        }
    }

    def B(x):
        return np.where(x < 660e3, cammarano_parameters[q_profile]["B"][0], cammarano_parameters[q_profile]["B"][1])

    def g(x):
        return np.where(x < 660e3, cammarano_parameters[q_profile]["g"][0], cammarano_parameters[q_profile]["g"][1])

    def a(x):
        return 0.2

    def omega(x):
        return 1.

    def Q_kappa(x):
        return np.where(x < 660e3, 1e3, 1e4)

    return gdrift.CammaranoAnelasticityModel(B=B, g=g, a=a, solidus=solidus, Q_bulk=Q_kappa, omega=omega)


def build_solidus():
    # Defining the solidus curve for manlte
    # First load the solidus curve of Andrault et al 2011 EPSL
    andrault_solidus = gdrift.RadialEarthModelFromFile(
        model_name="1d_solidus_Andrault_et_al_2011_EPSL",
        description="Andrault et al. 2011, EPSL")

    # Next load the solidus curve of Hirschmann 2000
    hirsch_solidus = gdrift.HirschmannSolidus()

    # Combining the two
    my_depths = []
    my_solidus = []

    for solidus_model in [hirsch_solidus, andrault_solidus]:
        # Getting minimum maximum of the profiles to re-discretise the profile
        d_min, d_max = solidus_model.min_max_depth("solidus temperature")
        dpths = np.arange(d_min, d_max, 10e3)

        # Add the values for our solidus curve
        my_depths.extend(dpths)
        my_solidus.extend(solidus_model.at_depth("solidus temperature", dpths))

    # Since we might have values outside the range of the solidus curve, we are better off with extrapolating
    ghelichkhan_et_al = SplineProfile(
        depth=np.asarray(my_depths),
        value=np.asarray(my_solidus),
        extrapolate=True,
        name="Ghelichkhan et al 2021")

    return ghelichkhan_et_al


def plot_layer(lons, lats, data, projection='Robinson', figsize=(12, 8),
               cmap='cmc.vik_r', title=None, colorbar_label='dV (%)',
               vmin=None, vmax=None, levels=20, save_path=None):
    """
    Plot a global map of seismic velocity data using Cartopy.

    Parameters:
    -----------
    lons : array-like
        Longitude values (degrees)
    lats : array-like
        Latitude values (degrees)
    data : array-like
        Data values to plot (e.g., velocity perturbations)
    projection : str, optional
        Cartopy projection name. Options: 'PlateCarree', 'Robinson',
        'Mollweide',
        'Orthographic', 'NorthPolarStereo', 'SouthPolarStereo'
    figsize : tuple, optional
        Figure size (width, height) in inches
    cmap : str, optional
        Colormap name
    title : str, optional
        Plot title
    colorbar_label : str, optional
        Label for the colorbar
    vmin, vmax : float, optional
        Color scale limits. If None, uses data min/max
    levels : int, optional
        Number of contour levels
    save_path : str or Path, optional
        Path to save the figure

    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """

    # Set up the projection
    projection_dict = {
        'PlateCarree': ccrs.PlateCarree(),
        'Robinson': ccrs.Robinson(),
        'Mollweide': ccrs.Mollweide(),
        'Orthographic': ccrs.Orthographic(central_longitude=0,
                                          central_latitude=0),
        'NorthPolarStereo': ccrs.NorthPolarStereo(),
        'SouthPolarStereo': ccrs.SouthPolarStereo()
    }

    proj = projection_dict.get(projection, ccrs.PlateCarree())

    # Create figure and axis
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    # Add map features
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3)
    ax.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.3)
    ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.3)

    # Set global extent for most projections
    if projection not in ['NorthPolarStereo', 'SouthPolarStereo']:
        ax.set_global()
        ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                     y_inline=False, linewidth=0.5, alpha=0.5)

    # Convert data to 2D arrays if needed
    lons = np.asarray(lons)
    lats = np.asarray(lats)
    data = np.asarray(data)

    # Handle different input formats
    if lons.ndim == 1 and lats.ndim == 1:
        # Create meshgrid if 1D arrays are provided
        lon_grid, lat_grid = np.meshgrid(lons, lats)
        if data.shape != lon_grid.shape:
            # Reshape data to match grid if necessary
            if data.size == lon_grid.size:
                data = data.reshape(lon_grid.shape)
            else:
                raise ValueError(
                    "Data dimensions don't match coordinate dimensions")
    else:
        # Use arrays as provided
        lon_grid, lat_grid = lons, lats

        # --- symmetric limits & levels about 0 ---
        A = float(np.nanmax(np.abs(data)))
        if not np.isfinite(A) or A == 0.0:
            A = 1.0  # tiny span for constant layers; keeps TwoSlopeNorm happy
        vmin, vmax = -A, +A

        # symmetric levels (e.g., 41 ticks from -A to +A)
        n_levels = 41  # tweak as you like; odd number keeps 0 exactly on a level
        levels = np.linspace(vmin, vmax, n_levels)

        # norm centered at 0; if something odd, fall back gracefully
        try:
            norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
        except ValueError:
            norm = Normalize(vmin=vmin, vmax=vmax)

    # Use contourf for filled contours
    contour = ax.contourf(lon_grid, lat_grid, data,
                          levels=levels, cmap=cmap, norm=norm,
                          transform=ccrs.PlateCarree(), extend='both')

    # Add colorbar
    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal',
                        pad=0.08, shrink=0.8, aspect=30)
    cbar.set_label(colorbar_label, fontsize=12)

    # Set title
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')

    # Adjust layout
    plt.tight_layout()

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")

    return fig, ax


#
if __name__ == "__main__":
    __main__(should_plot_layer=True)
