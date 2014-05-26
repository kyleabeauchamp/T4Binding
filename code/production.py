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

pdb_filename = "./equil_npt/%s_%s_%s.pdb" % (code, ff_name, water_name)
dcd_filename = "./production_hot/%s_%s_%s.dcd" % (code, ff_name, water_name)
log_filename = "./production_hot/%s_%s_%s.log" % (code, ff_name, water_name)

padding = 0.9 * u.nanometers
cutoff = 0.95 * u.nanometers
output_frequency = 5000
n_steps = 5000000

temperature = 300. * u.kelvin
pressure = 1.0 * u.atmospheres

pdb = app.PDBFile(pdb_filename)
ff = app.ForceField(which_forcefield, "benzene.xml", which_water)

topology = pdb.topology
positions = pdb.positions
system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)


desired_temperature = 750.
hot_atoms = np.arange(2638, 2650) - 1
rest.REST.perturb_system(system, temperature=desired_temperature, reference_temperature=temperature, hot_atoms=hot_atoms)

integrator = mm.LangevinIntegrator(temperature, 0.25 / u.picoseconds, 2.0 * u.femtoseconds)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))

simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature)


print("Using platform %s" % simulation.context.getPlatform().getName())
simulation.reporters.append(app.DCDReporter(dcd_filename, output_frequency))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True))
simulation.step(n_steps)
