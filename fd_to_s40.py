""" Adopted from the script Tom New has provided
"""
# %%
import numpy as np
from scipy.spatial import KDTree
import pyvista as pv
from pathlib import Path
import gdrift
import spherical_tools as st


pvtu_file = Path("/Users/sghelichkhani/Data/ADJOINT_2025/output/output_0.pvtu")

depth_res = 50  # Distance between depth slices in km
lon_res, lat_res = (
    360 // 1 + 1,
    180 // 1 + 1,
)  # Grid resolution (number of points in lon/lat)
cmb_depth = 2890  # CMB depth in km
r_earth_km = 6370  # Radius of the Earth in km

# path to models directory in S40RTS code
name = "HTZ22"
output_dir = Path(f"./interpolated_layers")
output_dir.mkdir(parents=True, exist_ok=True)


def load_pvtu_points(pvtu_file):
    model = pv.read(pvtu_file).clean()
    model.points *= r_earth_km * 1.e3 / 2.208  # normalise the model
    # drop unneeded arrays
    for array_name in model.point_data.keys():
        if array_name not in ["FullTemperature_CG", "Temperature_Deviation_CG"]:
            del model.point_data[array_name]
    # calculate T and T_av, dropping arrays after they become unneeded
    model.point_data["T"] = model["FullTemperature_CG"] * 3700 + 300
    model.point_data["dT"] = model["Temperature_Deviation_CG"] * (
        np.max(model["T"]) - np.min(model["T"])
    )
    model.point_data["T_av"] = model["T"] - model["dT"]
    # depth in meters
    model.point_data["depth"] = r_earth_km * 1.e3 - np.linalg.norm(
        model.points, axis=1
    )

    # initialise thermodynamic model
    slb_pyrolite = gdrift.ThermodynamicModel(
        "SLB_16",
        "pyrolite",
        temps=np.linspace(300, 4000),
        depths=np.linspace(0, cmb_depth * 1.e3),
    )

    # A temperautre profile representing the mantle average temperature
    # This is used to anchor the regularised thermodynamic table (we make sure the seismic speeds are the same at those temperature for the regularised and unregularised table)
    temperature_spline = gdrift.SplineProfile(
        depth=np.asarray([0., 500.e3, 2700.e3, 3000.e3]),
        value=np.asarray([300., 1000., 3000., 4000.]),
    )

    # Regularising the table
    # Regularisation works by saturating the minimum and maximum of variable gradients with respect to temperature.
    # Default values are between -inf and 0.0; which essentialy prohibits phase jumps that would otherwise render
    # v_s/v_p/rho versus temperature non-unique.
    linear_slb_pyrolite = gdrift.mineralogy.regularise_thermodynamic_table(
        slb_pyrolite,
        temperature_spline,
        regular_range={"v_s": [-0.5, 0.0], "v_p": [-0.5, 0.0], "rho": [-0.5, 0.0]},
    )

    cammarano_q_model = "Q6"  # choose model from cammarano et al., 2003
    anelasticity = gdrift.CammaranoAnelasticityModel.from_q_profile(
        cammarano_q_model
    )  # Instantiate the anelasticity model
    # apply anelastic correction
    linear_anelastic_slb_pyrolite = gdrift.apply_anelastic_correction(
        linear_slb_pyrolite, anelasticity
    )

    # compute seismic velocities in P and S
    model.point_data["Vs"] = linear_anelastic_slb_pyrolite.temperature_to_vs(
        temperature=np.array(model["T"]), depth=np.array(model["depth"])
    )
    model.point_data["Vp"] = linear_anelastic_slb_pyrolite.temperature_to_vp(
        temperature=np.array(model["T"]), depth=np.array(model["depth"])
    )

    # compute layer average seismic velocities in P and S
    model.point_data["Vs_av"] = linear_slb_pyrolite.temperature_to_vs(
        temperature=np.array(model["T_av"]), depth=np.array(model["depth"])
    )
    model.point_data["Vp_av"] = linear_slb_pyrolite.temperature_to_vp(
        temperature=np.array(model["T_av"]), depth=np.array(model["depth"])
    )

    # compute fractional perturbations in P and S
    model.point_data["dVs"] = 100 * (model["Vs"] - model["Vs_av"]) / model["Vs_av"]
    model.point_data["dVp"] = 100 * (model["Vp"] - model["Vp_av"]) / model["Vp_av"]

    return model.points / 1.e3, model["dVs"], model["dVp"], model["depth"] / 1.e3


# --- SLICE GENERATION ---
def get_depths_by_number(N):
    depths = np.linspace(0, cmb_depth, N).tolist()
    return depths


def get_depths_by_resolution(depth_res):
    # calculate number of depth slices starting from 0 such that the first N-1 slices are `depth_res` apart and the last slice is at the CMB depth
    N_slices = int(np.ceil(cmb_depth / depth_res)) + 1
    # calculate depth slice values
    depths = [depth_res * i for i in range(N_slices - 1)] + [cmb_depth]
    return depths


# --- GRID GENERATION ---
def get_lonlat_grid(lon_res, lat_res):
    lons = np.linspace(-180, 180, lon_res)
    lats = np.linspace(-90, 90, lat_res)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    return lon_grid, lat_grid


# --- INTERPOLATION ---
def idw_interpolate(tree, values, query_points, k=1000, eps=1e-12):
    dists, idxs = tree.query(query_points, k=k)
    weights = 1.0 / (dists + eps)
    weights /= weights.sum(axis=1, keepdims=True)
    interpolated_values = []
    for val in values:
        interpolated_values.append(np.sum(val[idxs] * weights, axis=1))
    return interpolated_values


fd_points, dVs, dVp, fd_depths = load_pvtu_points(pvtu_file)

## Get array of Firedrake layer average depths
# reshape such that the first index is a layer of the mesh and the second index is the depths at each point in that layer
firedrake_depths = fd_depths.reshape(-1,129).T
# take the mean of the depths at each layer
firedrake_depths = np.mean(firedrake_depths, axis=1)
# flip the array to have the shallowest layer first
firedrake_depths = np.flip(firedrake_depths)

## Create points for interpolation
# create array of interpolation depths
interpolation_depths = get_depths_by_resolution(depth_res)
# create arrays for lon and lat interpolation locations
interpolation_lons, interpolation_lats = get_lonlat_grid(lon_res, lat_res)
# convert to Cartesian coordinates

# create array to hold interpolated values
interpolated_dVs = np.empty((len(interpolation_depths) - 1, interpolation_lons.size))
interpolated_dVp = np.empty((len(interpolation_depths) - 1, interpolation_lons.size))

for i in range(len(interpolation_depths) - 1):
    print(f"Processing depth slice {interpolation_depths[i]} to {interpolation_depths[i + 1]} km")
    ## Select points between depth boundaries (slices will be at the midpoint between depths)
    ## Create mask for depth slice
    # we need to ensure that we have at least 1 Firedrake layer between depths[i-1] and depths[i], and if not, widen it gradually
    top = interpolation_depths[i]
    bottom = interpolation_depths[i + 1]
    step = (bottom - top) / 4 # amount by which to expand the depth slice above and below the slice
    found = False
    mask = np.empty_like(fd_depths, dtype=bool)
    while not found:
        idx_top = np.searchsorted(firedrake_depths, top, side='right')
        idx_bottom = np.searchsorted(firedrake_depths, bottom, side='left')
        if idx_top < idx_bottom:
            # we have at least one Firedrake layer in the slice
            found = True
            # create mask
            mask = (fd_depths >= top) & (fd_depths < bottom)
        else:
            # we need to widen the depth slice
            print(f"No Firedrake layers found in [{top},{bottom}] km, expanding slice...")
            top -= step
            bottom += step
            print(f"...expanded depth slice to [{top},{bottom}] km")

    slice_points = fd_points[mask]
    slice_dVs = dVs[mask]
    slice_dVp = dVp[mask]
    ## Build KDTree and interpolate
    tree = KDTree(slice_points)
    # convert grid points to Cartesian coordinates
    query_points = np.column_stack(
        [
            (r_earth_km - interpolation_depths[i]) * np.ones(interpolation_lons.size),
            interpolation_lons.reshape(-1),
            interpolation_lats.reshape(-1),
        ]
    )
    query_points = st.geo2cart(query_points, degrees=True)

    # interpolate using IDW and add to carrying arrays
    interpolated_dVs[i], interpolated_dVp[i] = idw_interpolate(
        tree, [slice_dVs, slice_dVp], query_points
    )
    # print(f"Max dVs: {interpolated_dVs[i].max():.2f}%; min dVs: {interpolated_dVs[i].min():.2f}%")


# %%
# save depths
np.savetxt(output_dir / "depth_layers.dat", interpolation_depths, fmt='%04d', newline='\n')

# save interpolated dVs values
interpolation_lons_to_save = interpolation_lons.reshape(-1)
interpolation_lats_to_save = interpolation_lats.reshape(-1)

lon_b, lat_b, dVs_b, dVp_b = np.broadcast_arrays(
    interpolation_lons_to_save[np.newaxis, :],
    interpolation_lats_to_save[np.newaxis, :],
    interpolated_dVs,
    interpolated_dVp,
)

dVs_to_save = np.stack([lon_b, lat_b, dVs_b], axis=-1)
dVp_to_save = np.stack([lon_b, lat_b, dVp_b], axis=-1)

for i in range(len(interpolation_depths) - 1):
    output_file_dVs = output_dir / f"{name}.dvs.layer.{i+1:03d}.dat"
    # output_file_dVp = output_dir / f"{name}.dvp.layer.{i+1:03d}.dat"
    np.savetxt(output_file_dVs, dVs_to_save[i], fmt=['%.2f', '%.2f', '%.8f'], delimiter=' ', newline='\n')
    # np.savetxt(output_file_dVp, dVp_to_save[i], fmt=['%.2f', '%.2f', '%.8f'], delimiter=' ', newline='\n')
