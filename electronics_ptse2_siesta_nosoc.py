""" IMPORTS """
import shelve
import os
from magic_dorye.actions.parser_actions import read_vasp_structure
from magic_dorye.actions.qe_actions.qe_creators import (create_qe_charge_density,create_qe_scf)
from magic_dorye.objects.material.atomic_structure import AtomicStructure
from magic_dorye.objects.material.band_structure import BandStructure
from magic_dorye.objects.material.density_of_states import DensityOfStates
from magic_dorye.objects.material.material import Material
from magic_dorye.objects.models import DFTModel
from magic_dorye.objects.computations.qe_calculation import QECalculation
from magic_dorye.objects.material.spin_physics import SpinPhysics
from magic_dorye.utils.file_manager import create_dir
from magic_dorye.objects.simulation import ExpDescription
from magic_dorye.actions.convergence_actions.scanning_tools import scan_parameter
from magic_dorye.actions.structure_actions.structure_transformations import crystal_to_cartesian, replicate_structure, remove_atom
from magic_dorye.plots.plot_output_basic import visualise_band_structure,visualise_1D_plot,visualise_dictionary
from magic_dorye.plots.plot_structure_advance import visualise_atomic_structure
from magic_dorye.actions.parser_actions.parser_qe_output import extract_qe_output_info
from magic_dorye.utils.pdos_tools import sum_qe_pdos
from magic_dorye.plots.plot_density import plot_density

""" USER DATA """
threads=4
np=12
structure_file = '/opt/project/dorye_projects/ptse2/ptse2_monolayer.vasp'

band_structure_data = {
    'kpath_label': 'hex'
}

spin_data = {
    'spin_polarisation': False,
    #'spin_configuration': [0.1, 0, 0],
    'spin_orbit_coupling': False,
    #'magnetic_species': ['Cr']
}

dft_data = {
    'pp_dir': '/opt/project/support/potentials/Pseudos/apeiron_pseudos',          
    'kpoints': [3,3, 1],   
    'wfc_cutoff': 50,
    'rho_cutoff': 500,
}

qe_data = {
'verbosity': 'high',
'outdir': './tmp',
'mixing_param': 0.1,
'mixing_mode': 'plain',
'occupations': 'smearing',
'dos_occupations': 'tetrahedra',
'smearing_type': 'cold',
'degauss': 0.015,
'force_calculation': False,
'pressure_calculation': False,
'conv_thr': 1e-5,
}

""" OBJECT GENERATION"""
structure_data = read_vasp_structure(structure_file)
atomic_structure = AtomicStructure.model_validate(structure_data)

dir_name= '/opt/project/dorye_projects/ptse2'
output = os.path.join(dir_name,'unit_cell.html')
visualise_atomic_structure(atomic_structure,output,1)
new_atomic_structure = replicate_structure(atomic_structure,(4,4,1))
output = os.path.join(dir_name,'super_cell.html')
visualise_atomic_structure(new_atomic_structure,output,1)

#Vacancy
dir_name= '/opt/project/dorye_projects/ptse2/qe/vacancies'
create_dir(dir_name)

new_atomic_structure = replicate_structure(atomic_structure,(2,2,1))
new_atomic_structure = remove_atom(new_atomic_structure,4)
new_atomic_structure = remove_atom(new_atomic_structure,5)
output = os.path.join(dir_name,'vacancy.html')
visualise_atomic_structure(new_atomic_structure,output,3)

dir_name= '/opt/project/dorye_projects/ptse2/qe/vacancies/charge_density'
create_dir(dir_name)

atomic_structure = AtomicStructure.model_validate(new_atomic_structure)
spin_physics = SpinPhysics.model_validate(spin_data)
band_structure = BandStructure.model_validate(band_structure_data)
material = Material(atomic_structure=atomic_structure,
                    spin_physics=spin_physics,
                    band_structure=band_structure)
dft_model = DFTModel.model_validate(dft_data) 
qe_calculation = QECalculation.model_validate(qe_data) 
exp_qe=ExpDescription(material=material,
                     model=dft_model,
                     calculation=qe_calculation,exp_name='ptse2_mono',threads=threads,np=np)




scf_file_name = os.path.join(dir_name, f'{exp_qe.exp_name}.scf.in')
exp_qe = create_qe_scf(scf_file_name,exp_qe)
exp_qe = create_qe_charge_density(scf_file_name,exp_qe)
qe_output_properties = extract_qe_output_info(exp_qe)
file_name = os.path.join(dir_name, f'{exp_qe.exp_name}.cd.xsf')
plot_density(file_name=file_name,isovalue=0.003, bond_threshold=2.0, rep_factors=(2,2,1))