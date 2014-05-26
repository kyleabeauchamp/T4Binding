import mdtraj as md
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
from repex import rest


code = "benzene"
ff_name = "amber99sbildn"
water_name = 'tip3p'

which_forcefield = "%s.xml" % ff_name
which_water = '%s.xml' % water_name

ligand_filename = "./benzene2.pdb"
out_pdb_filename = "./box/%s_%s_%s.pdb" % (code, ff_name, water_name)
dcd_filename = "./box/%s_%s_%s.dcd" % (code, ff_name, water_name)
log_filename = "./box/%s_%s_%s.log" % (code, ff_name, water_name)

padding = 0.9 * u.nanometers
cutoff = 0.95 * u.nanometers
output_frequency = 1000

n_equil_steps = 15000
n_steps = 500000


ligand_traj = md.load(ligand_filename)
ligand_xyz = ligand_traj.openmm_positions(0)
ligand_top = ligand_traj.top.to_openmm()

ff = app.ForceField(which_forcefield, "benzene.xml", which_water)

temperature = 300.
pressure = 1.0 * u.atmospheres

modeller = app.modeller.Modeller(ligand_top, ligand_xyz)
modeller.addSolvent(ff, padding=padding, model='tip3p')

topology = modeller.topology
positions = modeller.positions

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

desired_temperature = 700. * u.kelvin
hot_atoms = np.arange(12)
rest.REST.perturb_system(system, temperature=desired_temperature, reference_temperature=temperature, hot_atoms=hot_atoms)

system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))

integrator = mm.LangevinIntegrator(temperature, 1.0 / u.picoseconds, 1.0 * u.femtoseconds)
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
print('Minimizing...')
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature)
print('Equilibrating...')
simulation.step(n_equil_steps)

print('Production...')
simulation.reporters.append(app.DCDReporter(dcd_filename, output_frequency))
simulation.reporters.append(app.PDBReporter(out_pdb_filename, n_steps - 1))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True))

simulation.step(n_steps)
