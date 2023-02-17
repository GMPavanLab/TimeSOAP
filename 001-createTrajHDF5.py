from MDAnalysis import Universe as mdaUniverse
from MDAnalysis import transformations
from MDAnalysis.tests.datafiles import TPR, XTC
from HDF5er import MDA2HDF5
from sys import argv
from os import path
import re

getT = re.compile("ice([0-9]*)")

def createArgs(filename: str):
    T = getT.search(filename).group(1)
    return dict(T=T)

def getName(filename: str):
    topo = getT.search(filename).group(1)
    return topo

def worker(filename: str, ff: str, ts: str, frames="1ps"):
    filePosition = filename
    filename = path.basename(filename)
    # name = getName(filename)
    trajname = filename.split(".")[0]
    topo = "ice0.gro"  #GRO file
    args = createArgs(filename)
    args["frames"] = frames
    args["ts"] = ts
    args["ForceField"] = ff
    print(
        f"from trajectory {filename} and topology {topo}, creating a new trajectory hdf5 in"
    )
    print(f'"{trajname}.hdf5/Trajectories/{trajname}"')
    print(f"parameters:", args)
    u = mdaUniverse(topo, filename)
    print(len(u.atoms))
    # u.atoms.types = ["OW", "HW1", "HW2"] * int(len(u.atoms)/4)
    h2o = u.select_atoms("type O or type H")
    print(h2o)
    MDA2HDF5(h2o, trajname + ".hdf5", trajname, trajChunkSize=1000, attrs=args)
    # ref=u.select_atoms("index 0:1")
    # ref = mdaUniverse(topo)
    # NO FITTING -- SOAP is already translational-rotational invariant
    # u.trajectory.add_transformations(transformations.fit_rot_trans(u, ref))
    # MDA2HDF5(u, trajname + "_fitted.hdf5", f"{trajname}_fitted", trajChunkSize=1000)



if __name__ == "__main__":
    dumps = argv[1:]
    ts = "5fs"
    ff = "TIP4P/ICE"
    for d in dumps:
        worker(d, ff, ts)