from MDAnalysis import Universe as mdaUniverse
from MDAnalysis import transformations
from MDAnalysis.tests.datafiles import TPR, XTC
from HDF5er import MDA2HDF5
from sys import argv
from os import path

filename = argv[1]  # pass TRJ file (.xtc, .lammpsdump)
topo = argv[2]  # pass topology file (.gro, .data)

filename = path.basename(filename)
trajname = filename.split(".")[0]
print(
    f"from trajectory {filename} and topology {topo}, creating a new trajectory hdf5 in"
)
print(f'"{trajname}.hdf5/Trajectories/{trajname}"')
u = mdaUniverse(topo, filename)
print(len(u.atoms))
h2o = u.select_atoms("type O or type H")  # for ICE/WATER system
print(h2o)
MDA2HDF5(h2o, trajname + ".hdf5", trajname, trajChunkSize=1000)
# ref=u.select_atoms("index 0:1")
# ref = mdaUniverse(topo)
# NO FITTING -- SOAP is already translational-rotational invariant
# u.trajectory.add_transformations(transformations.fit_rot_trans(u, ref))
# MDA2HDF5(u, trajname + "_fitted.hdf5", f"{trajname}_fitted", trajChunkSize=1000)
