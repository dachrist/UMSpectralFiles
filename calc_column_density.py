# Python routine to calculate the column density of a chemical species using an ATMO
# chemistry output file
import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp1d
import sys

__all__ = ["get_column_density"]

# Reference and surface pressures - in Pa
p_ref = 1e4
p_surf = 2e7

# Molecular masses
molmass = {
    "H": 1.00784,
    "H2": 2.01588,
    "H2O": 18.01528,
    "NH3": 17.03052,
    "CH4": 16.0425,
    "CO": 28.0101,
    "CO2": 44.0095,
    "He": 4.0026020,
    "Na": 22.989769280,
    "K": 39.09830,
    "Li": 6.9410,
    "Cs": 132.90545190,
    "Rb": 85.46780,
    "TiO": 63.8664,
    "VO": 66.94090,
    "HCN": 27.0253,
}

#
def read_atmo_chem(fname):

    # Read atmo chemistry file
    data = nc.Dataset(fname)

    # Get pressure, abundances and molecule names
    pressure = data.variables["pressure"][:] * 0.1  # Convert to Pa
    abundances = data.variables["abundances"][:, :]
    names = data.variables["molname"][:]

    return pressure, abundances, names


# Function to calculate the column density
def calc_column_density_atmo(
    sp, abundance, gravity, pressure, abundance_h2, abundance_he
):

    # Get abundances of requested molecule, and H2 and He, at p_ref
    fint = interp1d(pressure, abundance)
    amol = fint(p_ref)

    fint = interp1d(pressure, abundance_h2)
    ah2 = fint(p_ref)

    fint = interp1d(pressure, abundance_he)
    ahe = fint(p_ref)

    # Calculate mean molecular weight (assume H2+He dominated)
    mean_mol_weight = ah2 * molmass["H2"] + ahe * molmass["He"]
    # Calculate column density
    sigma = (p_surf / gravity) * amol * molmass[sp] / mean_mol_weight

    return sigma


# Function to return the index of a species in the ATMO chemistry file
def imol(mol, mols):

    index = None

    for i in range(len(mols)):
        if mol == mols[i]:
            index = i

    if index is None:
        raise Exception("Error: molecules not included in ATMO file: ", mol)

    return index


# Main function
def get_column_density(
    fname=None,
    gravity=None,
    species=[
        "H2",
        "He",
        "H",
        "H2O",
        "CH4",
        "CO",
        "CO2",
        "NH3",
        "Na",
        "K",
        "Li",
        "Rb",
        "Cs",
        "TiO",
        "VO",
        "HCN",
    ],
):

    col_out = {}
    # Check required inputs are defined
    if fname is None:
        raise Exception("Error: fname is not defined")

    if gravity is None:
        raise Exception("Error: you must define gravity in ms-2")

    # Read chemistry file
    pressure, abundances, names = read_atmo_chem(fname)

    # Get indices of requested species, clean up formatting
    species_name = []
    for j in range(len(names[:, 0])):
        name = (
            names[j, 0]
            + names[j, 1]
            + names[j, 2]
            + names[j, 3]
            + names[j, 4]
            + names[j, 5]
            + names[j, 6]
            + names[j, 7]
            + names[j, 8]
            + names[j, 9]
        )
        species_name.append((name.decode("utf-8")).strip())

    abundance = np.zeros(len(pressure))
    abundance_h2 = abundances[imol("H2", species_name), :]
    abundance_he = abundances[imol("He", species_name), :]

    # Loop over molecules, calculate and print column density
    for sp in species:
        molmass_sp = molmass[sp]
        abundance = abundances[imol(sp, species_name), :]
        sigma = calc_column_density_atmo(
            sp, abundance, gravity, pressure, abundance_h2, abundance_he
        )
        col_out[sp] = sigma
    return col_out
