from SOAPify import saponifyGroup
import h5py
from sys import argv


def worker(trajFileName: str) -> None:
    soapFileName = trajFileName.split(".")[0] + "soap.hdf5"
    print(trajFileName, soapFileName)
    with h5py.File(trajFileName, "r") as workFile, h5py.File(
        soapFileName, "a"
    ) as soapFile:
        saponifyGroup(
            trajContainers=workFile["Trajectories"],
            SOAPoutContainers=soapFile.require_group("SOAP"),
            SOAPOutputChunkDim=1000,
            SOAPatomMask="O",  # to choose OW as SOAP center
            SOAPnJobs=32,
            SOAPrcut=10,
            SOAPnmax=8,
            SOAPlmax=8,
        )


if __name__ == "__main__":
    if len(argv) != 2:
        print(f"Usage: python3 {argv[0]} <trajHDF5file>")
        exit(1)
    worker(argv[1])
