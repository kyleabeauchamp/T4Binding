import statsmodels.api as sm
import mdtraj as md
import simtk.unit as u

code = "benzene"
ff_name = "amber99sbildn"
water_name = 'tip3p'

which_forcefield = "%s.xml" % ff_name
which_water = '%s.xml' % water_name

out_pdb_filename = "./water/box.pdb"


temperature = 1600.0 * u.kelvin
dcd_filename = "./water/%s_%s_%s_%s.dcd" % (code, ff_name, water_name, temperature)

t1 = md.load(dcd_filename, top=out_pdb_filename)
np.diff(t1.xyz[:, 0], 0).std(0)

temperature = 300.0 * u.kelvin
dcd_filename = "./water/%s_%s_%s_%s.dcd" % (code, ff_name, water_name, temperature)

t0 = md.load(dcd_filename, top=out_pdb_filename)
np.diff(t0.xyz[:, 0], 0).std(0)

d0 = md.compute_distances(t0, np.array([[0, 15]]))[:, 0]
d1 = md.compute_distances(t1, np.array([[0, 15]]))[:, 0]

plot(d0, label="300")
plot(d1, label="hot")
