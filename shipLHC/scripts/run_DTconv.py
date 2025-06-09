#!/usr/bin/env python
import os
import resource
from argparse import ArgumentParser

import ROOT
import shipRoot_conf

# for the geometry
import SndlhcGeo


def mem_monitor():
    # Getting virtual memory size
    pid = os.getpid()
    with open(os.path.join("/proc", str(pid), "status")) as f:
        lines = f.readlines()
    _vmsize = [line_i for line_i in lines if line_i.startswith("VmSize")][0]
    vmsize = int(_vmsize.split()[1])
    # Getting physical memory size
    pmsize = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(
        "memory: virtuell = %5.2F MB  physical = %5.2F MB"
        % (vmsize / 1.0e3, pmsize / 1.0e3)
    )


firstEvent = 0

shipRoot_conf.configure()

parser = ArgumentParser()
parser.add_argument(
    "-f", "--inputFile", dest="inputFile", help="single input file", required=True
)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True) # maybe to extract some params?
parser.add_argument(
    "-n",
    "--nEvents",
    dest="nEvents",
    type=int,
    help="number of events to process",
    default=1000000,
)
parser.add_argument("-d", "--Debug", dest="debug", help="debug", default=False)

options = parser.parse_args()
# -----Timer-------------
timer = ROOT.TStopwatch()
timer.Start()

# outfile name
tmp = options.inputFile.split("/")
# outFile = tmp[len(tmp) - 1].replace(".root", "_dig.root")
tmpOutFile = tmp[len(tmp) - 1].replace(".root", "_dig.root")
runN = options.inputFile.split("/")[-2][4:]
outFile = "/eos/user/l/lmozzina/MiniDT/matching/run_" + runN + "/" + tmpOutFile

# -----Create geometry----------------------------------------------
snd_geo = SndlhcGeo.GeoInterface(options.geoFile)

# if needed to read the etector geometry
lsOfGlobals  = ROOT.gROOT.GetListOfGlobals()
# DriftTubeDet     = lsOfGlobals.FindObject('DriftTube')

run = ROOT.FairRunAna()
ioman = ROOT.FairRootManager.Instance()
# ioman.RegisterInputObject("DriftTube", snd_geo.modules["DriftTube"])
# Set input
fileSource = ROOT.FairFileSource(options.inputFile)
run.SetSource(fileSource)
# Set output
outFile_sink = ROOT.FairRootFileSink(outFile)
run.SetSink(outFile_sink)

# Set number of events to process
inRootFile = ROOT.TFile.Open(options.inputFile)
# # raw SND
# inTree = inRootFile.Get("data")  # FIXME input tree name
# converted SND
inTree = inRootFile.Get("rawConv")  # FIXME input tree name
nEventsInFile = inTree.GetEntries()
nEvents = min(nEventsInFile, options.nEvents)

rtdb = run.GetRuntimeDb()

# # raw SND
# ioman.RegisterInputObject("data", inRootFile)
# converted SND
ioman.RegisterInputObject("rawConv", inRootFile)

run.SetEventHeaderPersistence(False)
xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()

ConvDriftTubeTask = ROOT.ConvDriftTubeRawData()
run.AddTask(ConvDriftTubeTask)
run.Init()
run.Run(firstEvent, nEvents)
ConvDriftTubeTask.PrintMatchedEntries()

# outTree = ioman.GetOutTree()
# treeList = ROOT.TList()
# # treeList.Add(inTree)
# treeList.Add(outTree)
# # mergedTree = ROOT.TTree("new_tree", "new_tree")
# mergedTree = inTree.CloneTree()
# mergedTree.Merge(treeList)
# mergedTree.Write()
# # mergedTree = ROOT.TTree.MergeTrees(treeList)
# print(mergedTree.Print())

timer.Stop()
rtime = timer.RealTime()
ctime = timer.CpuTime()
print(" ")
print("Real time ", rtime, " s, CPU time ", ctime, "s")
