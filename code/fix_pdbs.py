import simtk.openmm.app as app
import pdbfixer

fixer = pdbfixer.PDBFixer(pdbid='3DMX')
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.removeHeterogens(True)
fixer.addMissingHydrogens(7.0)
numChains = len(list(fixer.topology.chains()))
fixer.removeChains(range(1, numChains))
app.PDBFile.writeFile(fixer.topology, fixer.positions, open("./pdb_fixed/3DMX.pdb", 'w'))


"""
 cat ~/pdb/3dmx.pdb|grep BNZ > native_benzene.pdb
reduce native_benzene.pdb  > out.pdb

"""

import mdtraj as md
import mdtraj.geometry.alignment

t0 = md.load("./benzene.pdb")
t0.xyz -= t0.xyz.mean(0).mean(0)
t0.restrict_atoms(arange(5))

t1 = md.load("./native_benzene.pdb")
t1.restrict_atoms(arange(6))
t1.xyz -= t1.xyz.mean(0).mean(0)
t1.restrict_atoms(arange(5))

trans = md.geometry.alignment.compute_transformation(t0.xyz[0], t1.xyz[0])

t0 = md.load("./benzene.pdb")
t0.xyz -= t0.xyz.mean(0).mean(0)
t0.xyz[0] = trans.transform(t0.xyz[0])


t1 = md.load("./native_benzene.pdb")

t0.xyz += t1.xyz.mean(0).mean(0)


t0.save("./benzene2.pdb")
