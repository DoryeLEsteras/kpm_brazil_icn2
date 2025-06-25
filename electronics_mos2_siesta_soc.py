"""IMPORTS"""

import shelve
import os
from magic_dorye.actions.parser_actions import read_vasp_structure
from magic_dorye.actions.siesta_actions.siesta_creators import (
    create_siesta_bands,
    create_siesta_scf,
)
from magic_dorye.objects.material.atomic_structure import AtomicStructure
from magic_dorye.objects.material.band_structure import BandStructure
from magic_dorye.objects.material.density_of_states import DensityOfStates
from magic_dorye.objects.material.optical_properties import OpticalProperties
from magic_dorye.objects.material.material import Material
from magic_dorye.objects.models import DFTModel
from magic_dorye.objects.computations.siesta_calculation import SiestaCalculation
from magic_dorye.objects.material.spin_physics import SpinPhysics
from magic_dorye.utils.file_meneger import create_dir
from magic_dorye.objects.exp_description import ExpDescription
from magic_dorye.actions.convergence_actions.scanning_tools import scan_parameter
from magic_dorye.actions.structure_actions.structure_transformations import (
    crystal_to_cartesian,
    replicate_structure,
    remove_atom,
)
from magic_dorye.plots.plot_output_basic import visualise_dictionary
from magic_dorye.plots.plot_structure_advance import visualise_atomic_structure
from magic_dorye.actions.parser_actions.parser_siesta_output import (
    extract_siesta_output_info,
)

""" USER DATA """
threads = 1
np = 10
structure_file = "mos2_monolayer.vasp"

band_structure_data = {"kpath_label": "hex"}

spin_data = {
    "spin_polarisation": False,
    "spin_orbit_coupling": True,
}

dft_data = {
    "pp_dir": "/opt/project/dorye_projects/brazil_kpm/lda_pp",
    "kpoints": [12, 12, 1],
    "dos_kpoints": [20, 10, 1],
    "mesh_cutoff": (1200, "Ry"),
    "orbital_basis": "DZP",
    "electronic_temperature": (50, "K"),
}

siesta_data = {
    "mixing_param": 0.1,
    "mixing_mode": "Broyden",
    "mixing_history": 10,
    "max_scf_iterations": 800,
}

""" OBJECT GENERATION"""
structure_data = read_vasp_structure(structure_file)

atomic_structure = AtomicStructure.model_validate(structure_data)
spin_physics = SpinPhysics.model_validate(spin_data)
band_structure = BandStructure.model_validate(band_structure_data)

material = Material(
    atomic_structure=atomic_structure,
    spin_physics=spin_physics,
    band_structure=band_structure,
)

dft_model = DFTModel(**dft_data)
siesta_calculation = SiestaCalculation(**siesta_data)

exp_siesta = ExpDescription(
    material=material,
    model=dft_model,
    calculation=siesta_calculation,
    exp_name="mos2",
    threads=threads,
    np=np,
)

dir_name = "/opt/project/dorye_projects/brazil_kpm/mos2/siesta/soc/bands"
create_dir(dir_name)
scf_file_name = os.path.join(dir_name, f"{exp_siesta.exp_name}.fdf")
exp_siesta = create_siesta_scf(scf_file_name, exp_siesta)
# exp_siesta = create_siesta_bands(scf_file_name, exp_siesta)
# siesta_output_properties = extract_siesta_output_info(exp_siesta)
# output = os.path.join(dir_name, "bands.html")
# visualise_band_structure(
#    exp_siesta.material.computed_info["Siesta"]["band_structure"],
#    output_file=output,
#    fermi=siesta_output_properties["fermi_energy"],
#    x_labels=["Γ", "M", "K", "Γ"],
#    high_sym_points=[],
# )
