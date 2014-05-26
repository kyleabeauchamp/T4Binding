import mdtraj as md
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

code = "3DMX"
ff_name = "amber99sbildn"
water_name = 'tip3p'

which_forcefield = "%s.xml" % ff_name
which_water = '%s.xml' % water_name

pdb_filename = "./pdb_fixed/%s.pdb" % code
ligand_filename = "./benzene2.pdb"
out_pdb_filename = "./equil_npt/%s_%s_%s.pdb" % (code, ff_name, water_name)
dcd_filename = "./equil_npt/%s_%s_%s.dcd" % (code, ff_name, water_name)
log_filename = "./equil_npt/%s_%s_%s.log" % (code, ff_name, water_name)

padding = 0.9 * u.nanometers
cutoff = 0.95 * u.nanometers
output_frequency = 1000
n_steps = 50000

trj = md.load(pdb_filename).stack(md.load(ligand_filename))

ff = app.ForceField(which_forcefield, "benzene.xml", which_water)

temperature = 300.
pressure = 1.0 * u.atmospheres


modeller = app.Modeller(trj.top.to_openmm(), trj.openmm_positions(0))
modeller.addSolvent(ff, padding=padding, model='tip3p')

topology = modeller.topology
positions = modeller.positions

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, 1.0 / u.picoseconds, 1.0 * u.femtoseconds)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(temperature)
print('Equilibrating...')

simulation.step(25000)  # Don't even save the first 25 ps

simulation.reporters.append(app.DCDReporter(dcd_filename, output_frequency))
simulation.reporters.append(app.PDBReporter(out_pdb_filename, n_steps - 1))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True))
simulation.step(n_steps)
