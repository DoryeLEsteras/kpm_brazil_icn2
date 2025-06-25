"""IMPORTS"""

import shelve
import os
from magic_dorye.actions.parser_actions import read_vasp_structure
from magic_dorye.actions.qe_actions.qe_creators import (
    create_qe_bands,
    create_qe_pdos,
    create_qe_scf,
)
from magic_dorye.objects.material.atomic_structure import AtomicStructure
from magic_dorye.objects.material.band_structure import BandStructure
from magic_dorye.objects.material.density_of_states import DensityOfStates
from magic_dorye.objects.material.material import Material
from magic_dorye.objects.models import DFTModel
from magic_dorye.objects.computations.qe_calculation import QECalculation
from magic_dorye.objects.material.spin_physics import SpinPhysics
from magic_dorye.utils.file_meneger import create_dir
from magic_dorye.objects.exp_description import ExpDescription
from magic_dorye.actions.convergence_actions.scanning_tools import scan_parameter
from magic_dorye.actions.structure_actions.structure_transformations import (
    crystal_to_cartesian,
    replicate_structure,
    remove_atom,
)
from magic_dorye.plots.plot_output_basic import (
    visualise_band_structure,
    visualise_1D_plot,
    visualise_dictionary,
)
from magic_dorye.plots.plot_structure_advance import visualise_atomic_structure
from magic_dorye.actions.parser_actions.parser_qe_output import extract_qe_output_info
from magic_dorye.utils.pdos_tools import sum_qe_pdos

""" USER DATA """
threads = 1
np = 20
structure_file = "mos2_monolayer.vasp"

band_structure_data = {"kpath_label": "hex"}

spin_data = {
    "spin_polarisation": False,
    "spin_orbit_coupling": True,
}

dft_data = {
    "pp_dir": "/opt/project/dorye_projects/brazil_kpm/lda_pp",
    "kpoints": [12, 12, 1],
    "dos_kpoints": [20, 20, 1],
    "wfc_cutoff": 120,
    "rho_cutoff": 1200,
}

qe_data = {
    "verbosity": "high",
    "outdir": "./tmp",
    "mixing_param": 0.1,
    "mixing_mode": "plain",
    "occupations": "smearing",
    "dos_occupations": "tetrahedra",
    "smearing_type": "cold",
    "degauss": 0.003,
    "conv_thr": 1e-7,
    "force_calculation": False,
    "pressure_calculation": False,
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
qe_calculation = QECalculation(**qe_data)

exp_qe = ExpDescription(
    material=material,
    model=dft_model,
    calculation=qe_calculation,
    exp_name="mos2",
    threads=threads,
    np=np,
)

"""
dir_name = "/opt/project/dorye_projects/brazil_kpm/mos2/structure"
create_dir(dir_name)
output = os.path.join(dir_name, "unit_cell.html")
visualise_atomic_structure(material.atomic_structure, output, 1)

new_atomic_structure = replicate_structure(atomic_structure, (3, 3, 2))
output = os.path.join(dir_name, "super_cell.html")
visualise_atomic_structure(new_atomic_structure, output, 1)
"""

dir_name = "/opt/project/dorye_projects/brazil_kpm/mos2/qe/soc/bands"
create_dir(dir_name)
scf_file_name = os.path.join(dir_name, f"{exp_qe.exp_name}.scf.in")
exp_qe = create_qe_scf(scf_file_name, exp_qe)
exp_qe = create_qe_bands(scf_file_name, exp_qe)
qe_output_properties = extract_qe_output_info(exp_qe)
output = os.path.join(dir_name, "bands.html")
visualise_band_structure(
    exp_qe.material.computed_info["QE"]["band_structure"],
    output_file=output,
    fermi=qe_output_properties["fermi_energy"],
    x_labels=["Γ", "M", "K", "Γ"],
    high_sym_points=[0.0000, 0.5774, 0.9104, 1.5764],
)


dir_name = "/opt/project/dorye_projects/brazil_kpm/mos2/qe/soc/dos"
create_dir(dir_name)
scf_file_name = os.path.join(dir_name, f"{exp_qe.exp_name}.scf.in")
exp_qe = create_qe_scf(scf_file_name, exp_qe)
exp_qe = create_qe_pdos(scf_file_name, exp_qe)
qe_output_properties = extract_qe_output_info(exp_qe)
(energies, total_dos, pdos) = sum_qe_pdos(exp_qe)
output = os.path.join(dir_name, "total_dos.html")
visualise_1D_plot(
    x_values=energies,
    y_values_list=total_dos,
    fermi=qe_output_properties["fermi_energy"],
    output_file=output,
)
output = os.path.join(dir_name, "pdos.html")
visualise_dictionary(
    x_values=energies,
    y_values_dict=pdos,
    fermi=qe_output_properties["fermi_energy"],
    output_file=output,
)

