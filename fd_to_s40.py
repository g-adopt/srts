""" Adopted from the script Tom New has provided
"""
# %%
import matplotlib.pyplot as plt
import numpy as np
from gdrift.profile import SplineProfile
import pyvista as pv
from pathlib import Path
import gdrift
from scipy.spatial import cKDTree
# import spherical_tools as st


def __main__():
    # Dimensional constants for conversion and model
    dimensional_constants = get_dimensional_constants()

    # Path to the pvtu file -- this is a gadopt model
    pvtu_file = Path(
        "/Users/sghelichkhani/Data/ADJOINT_2025/output/output_0.pvtu")

    # path to models directory in S40RTS code
    name = "0Ma_C52_8e8_Kappa"
    output_dir = Path(f"./interpolated_layers_{name}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # load the model
    model = load_pvtu_points(pvtu_file)
    model, avg_temperature_spline = generate_thermodynamic_model(model)

    # From this point on we are interpolating into a depth profile consitent with
    # srts example model
    # create array of interpolation depths
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
    ).query(np.column_stack(((depth_grid).flatten()/1e3, lon_grid.flatten(), lat_grid.flatten())), k=100)

    # Compute Gaussian weights based on distances
    sigma = 1.0  # Gaussian width parameter - adjust as needed
    weights = np.exp(-0.5 * (dists / sigma)**2)
    weights /= np.sum(weights, axis=1, keepdims=True)  # normalize weights

    # Compute weighted average of dVs using einsum
    rts_dVs = np.einsum('ij,ij->i', weights,
                        model.point_data["dVs"][inds]).reshape(depth_grid.shape)

    # # save depths
    np.savetxt(output_dir / "depth_layers.dat",
               depth_profile/1e3, fmt='%04d', newline='\n')

    for ir, dpth in enumerate(depth_profile):
        output_file_dVs = output_dir / f"{name}.dvs.layer.{ir:03d}.dat"

        # output_file_dVp = output_dir / f"{name}.dvp.layer.{i+1:03d}.dat"
        with open(output_file_dVs, mode="w") as f:
            f.write("\n".join([f"{lat:+.2f} {lon:+.2f} {0 if ir in [0, len(depth_profile)-1] else dat:.5f}" for lat, lon, dat in zip(
                lat_grid[ir, :, :].flatten(), lon_grid[ir, :, :].flatten(), rts_dVs[ir, :, :].flatten())]))


def get_dimensional_constants():
    return {
        "T_0": 3700,
        "T_1": 300,
        "r_max": 2.208,
        "r_min": 1.208,
        "nlayers": 65,
        "n_lons": 361,
        "n_lats": 181,
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


def generate_thermodynamic_model(model):
    dimensional_constants = get_dimensional_constants()
    # This is the depth profile for which we generate the thermodynamic model
    depth_profile = np.linspace(
        0, gdrift.R_earth - gdrift.R_cmb, dimensional_constants["nlayers"])
    # initialise thermodynamic model
    slb_pyrolite = gdrift.ThermodynamicModel(
        "SLB_16", "pyrolite", depths=depth_profile,
    )

    # Computing a radial profile of average temperature
    dists, inds = cKDTree(model.point_data["depth"][:, np.newaxis]).query(
        depth_profile[:, np.newaxis], k=100)

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
        temperature=np.array(model["T"]), depth=np.array(model["depth"]))
    model.point_data["Vp"] = linear_anelastic_slb_pyrolite.temperature_to_vp(
        temperature=np.array(model["T"]), depth=np.array(model["depth"]))

    # compute layer average seismic velocities in P and S
    model.point_data["Vs_av"] = linear_slb_pyrolite.temperature_to_vs(
        temperature=np.array(model["T_av"]), depth=np.array(model["depth"]))
    model.point_data["Vp_av"] = linear_slb_pyrolite.temperature_to_vp(
        temperature=np.array(model["T_av"]), depth=np.array(model["depth"]))

    # compute fractional perturbations in P and S
    model.point_data["dVs"] = (
        100 * (model["Vs"] - model["Vs_av"]) / model["Vs_av"])
    model.point_data["dVp"] = (
        100 * (model["Vp"] - model["Vp_av"]) / model["Vp_av"])

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


#
if __name__ == "__main__":
    __main__()
