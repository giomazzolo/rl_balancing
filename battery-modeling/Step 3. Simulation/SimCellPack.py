import subprocess
import time
import mmapTx
import numpy as np

class SimCellPack:

    def __init__(self, name = None, format_type = "d", in_size = 1, out_size = 1, blocking = True):
        return

        
def getSimData(numCells):

    mmaptx = mmapTx.Mmaptx(name="simCell", format_type="d", in_size=numCells+1, out_size=numCells+1, blocking=True)

    data_list = []
    cmd = 1
    feedback = [0]*numCells # Creates array of size numCells initialized with 0s
    start = time.time()

    # Start recieving data
    while True:

        mmaptx.write(cmd, *feedback)
        data_in = mmaptx.read()
        
        if data_in[0] == False:
            break
        else:
            data_list.append(list(data_in[1:]))

        # Toy balancing feedback
        minSoc = min(data_in[1:])
        feedback = [x - minSoc for  x in data_in[1:]]
            

    # Send 
    cmd = 0 # Send stop command
    feedback = [0]*numCells
    mmaptx.write(cmd, *feedback)

    # Print out time
    end = time.time()
    print("Execution time:", end-start)

    # Clean and close mmap
    mmaptx.cleanup()
    del mmaptx

    return data_list

def runSim(cmd):
    # stdout=subprocess.PIPE, creationflags=0x08000000
    return subprocess.Popen(cmd, stdout=subprocess.PIPE)

def startSim(numCells, cycles, sampleFactor):

    cmd = "runPackSim " + str(numCells) + " " + str(cycles) + " \"..\\data\\us60Power.mat\" \"..\\data\\P14model_dynamic.mat\" [0,1,1,1,1,1] " + str(sampleFactor)

    process = runSim(cmd)
    data = getSimData(numCells)
    process.communicate()
    return data