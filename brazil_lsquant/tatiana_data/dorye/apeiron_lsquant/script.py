from magic_dorye.actions.scanning_actions.scanning_tools import scan_parameter
from magic_dorye.actions.lsquant_actions.lsquant_creators import create_lsquant_dos
from magic_dorye.objects.models.lsquant_model import LSQuantModel
from magic_dorye.objects.simulation import Simulation

np=1
threads=1
data = {
'supercell_dims': [400,400,1],
'broadening':0.05,
'n_moments':100,
'bounds':(-19,11),
}

lsquant_model = LSQuantModel(**data)

sim_lsquant = Simulation(
    model=lsquant_model,
    sim_name="mos2.titan",
    threads=threads,
    np=np,
)

dir_name = "/opt/project/dorye_projects/kpm_brazil_icn2/brazil_lsquant/tatiana_data/dorye/apeiron_lsquant"

create_lsquant_dos(dir_name,sim_lsquant)

#bounds = [(-100,100),(-90,90),(-80,80),(-70,70),(-60,60)]
#for bound in bounds:
#    print(bound)
#    sim_lsquant.model.bounds = bound
#    create_lsquant_dos(dir_name,sim_lsquant)
#simulations = scan_parameter('/opt/project/dorye_projects/kpm_brazil_icn2/brazil_lsquant/tatiana_data/dorye/apeiron_lsquant/BOUNDS',sim_lsquant,'bounds',
#                                 [create_lsquant_dos],bounds)