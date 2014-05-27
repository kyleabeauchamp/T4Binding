import mdtraj as md
import numpy as np
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
from repex import rest

code = "3DMX"
ff_name = "amber99sbildn"
water_name = 'tip3p'

which_forcefield = "%s.xml" % ff_name
which_water = '%s.xml' % water_name

pdb_filename = "./pdb_fixed/%s.pdb" % code
ligand_filename = "./benzene2.pdb"

padding = 0.9 * u.nanometers
cutoff = 0.95 * u.nanometers
output_frequency = 5000
n_steps = 50000000

temperature = 300. * u.kelvin
pressure = 1.0 * u.atmospheres

pdb = app.PDBFile(pdb_filename)
ff = app.ForceField(which_forcefield, "benzene.xml", which_water)

protein_traj = md.load(pdb_filename)
protein_top = protein_traj.top.to_openmm()
protein_xyz = protein_traj.openmm_positions(0)

ligand_traj = md.load(ligand_filename)
ligand_xyz = ligand_traj.openmm_positions(0)
ligand_top = ligand_traj.top.to_openmm()

modeller = app.modeller.Modeller(protein_top, protein_xyz)
modeller.add(ligand_top, ligand_xyz)
modeller.addSolvent(ff, padding=padding, model='tip3p')

topology = modeller.topology
positions = modeller.positions

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))

temperature_tuple = (500., 1600.)

hot_atom_lists = []
hot_atom_lists.append(np.arange(2637))
hot_atom_lists.append(np.arange(2638, 2650) - 1)

rest.REST.soften_system(system, temperature_tuple=temperature_tuple, reference_temperature=temperature, hot_atom_lists=hot_atom_lists)


integrator = mm.LangevinIntegrator(temperature, 1.0 / u.picoseconds, 1.0 * u.femtoseconds)
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
print('Minimizing...')
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature)
print('Production...')

dcd_filename = "./production/multi.dcd"
log_filename = "./production/multi.log"

print("Using platform %s" % simulation.context.getPlatform().getName())
simulation.reporters.append(app.DCDReporter(dcd_filename, output_frequency))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True))
simulation.step(n_steps)
