import mdtraj as md

code = "benzene"
ff_name = "amber99sbildn"
water_name = 'tip3p'

which_forcefield = "%s.xml" % ff_name
which_water = '%s.xml' % water_name

out_pdb_filename = "./box/%s_%s_%s.pdb" % (code, ff_name, water_name)
dcd_filename = "./water_hot/%s_%s_%s.dcd" % (code, ff_name, water_name)

t1 = md.load(dcd_filename, top=out_pdb_filename)
np.diff(t1.xyz[:, 0], 0).std(0)


dcd_filename = "./water/%s_%s_%s.dcd" % (code, ff_name, water_name)

t0 = md.load(dcd_filename, top=out_pdb_filename)
np.diff(t0.xyz[:, 0], 0).std(0)

d0 = md.compute_distances(t0, np.array([[0, 15]]))[:, 0]
d1 = md.compute_distances(t1, np.array([[0, 15]]))[:, 0]
