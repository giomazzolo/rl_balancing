import subprocess
import mmapTx
import os
import signal

import sys,os

cellModels = ["E1", "E2", "P14"]
        
balancingTypes = ["active", "passive"]

drive_cycle_profiles = ["sc03", "bcdc", "ftp", "hwfet", "nycc", "ucds", "us06", "udds"]

getSOCsWhenList = ["all", "charge", "discharge", "dis/charge", "rest"]

class SimCellPack:
# This class calls a cell simulation executable with the specified parameters
    # Must be called from within rl_balance
    
    def __init__(self, cellModel, numCells, simCycles, seed, profile, balancing, getSOCsWhen, sampleFactor, utilization):
        # Initialize the Cell Simulation class
        
        self.simSubProcess = None
        self.mmaptx = None

        self.numCells = numCells
        self.simCycles = simCycles

        self.seed = seed

        # Cell models available for the simulation
        # ["A12", "ATL", "SAM"] skipped cell models

        try:
            index = os.getcwd().split("\\").index("rl_balancing")
        except:
            print(("rl_balance folder not found, make sure to run this script inside rl_balancing folder."))
            exit

        self.modelsDir = "\\".join(os.getcwd().split("\\")[:index+1]) + "\\cell_simulation"
        
        self.simExecutable = "\"" + self.modelsDir + "\\runPackSim.exe\""

        if cellModel not in cellModels:
            print("Cell model specified not in the list: " + str(cellModels))
            exit()
        else:
            self.cellModelPath = "\"" + self.modelsDir + "\\data\\cell_models\\" + cellModel + "model.mat\""
        
        if profile not in drive_cycle_profiles:
            print("Drive cycle profile specified not in the list: " + str(drive_cycle_profiles))
            exit()
        else:        
            self.profilePath = "\"" + self.modelsDir + "\\data\\drive_cycle_profiles\\" + cellModel + "_" + profile + "_profile.mat\""

        if balancing not in balancingTypes:
            print("Balancing type specified not in the list: " + balancingTypes)
            exit()
        else:
            self.balancing = balancing

        if getSOCsWhen not in getSOCsWhenList:
            print("getSOCS type specified not in the list: " + getSOCsWhenList)
            exit()
        else:
            self.getSOCsWhen = getSOCsWhen

        self.cellRandOpts = "[0,1,1,1,1,1]"

        self.sampleFactor = sampleFactor

        self.utilization =  "\"" + str(utilization).replace(" ", "") + "\""
        
        
    def startSim(self):
    # Starts the execution of the cell executable simulator based on the configured parameters
        self.simCmd = 0
        
        self.cmdExe = " ".join([self.simExecutable, str(self.numCells), str(self.simCycles),        \
                                self.profilePath, self.cellModelPath, str(self.seed), self.balancing,    \
                                self.getSOCsWhen, self.cellRandOpts, str(self.sampleFactor), self.utilization])
                   
        self.simSubProcess = subprocess.Popen(self.cmdExe, stdout=subprocess.PIPE, shell=True)

        self.mmaptx = mmapTx.Mmaptx(name="simCell", format_type="d", in_size=self.numCells+1, out_size=self.numCells+1, blocking=True)


    def getSimStep(self):
    # Returns a simulation step data
        # The cell simulator will await for feedback before continuing the next step
        # The returns the simulation state and the data
            # If the simulation state is 1, that means that the siulation is still ongoing
            # If the simulation state is 0, that means that the simulation finished

        data_in = self.mmaptx.read()
        return data_in[0], data_in[1:]
    
    
    def sendSimFeedback(self, feedback = None):
    # Sends control feedback to the cell somulation
        self.mmaptx.write(self.simCmd, *feedback)

    def resetSim(self):
        self.simCmd = 2 # Reset code for the cell simulator
        self.mmaptx.write(self.simCmd) # Send reset command
        self.simCmd = 0
    
    def stopSim(self):
    # Stops the simulation
        if self.mmaptx != None:
            self.simCmd = 1.0
            self.mmaptx.write(self.simCmd)
            if self.simSubProcess != None:
                stdout, stderr = self.simSubProcess.communicate()
                try:
                    os.killpg(os.getpgid(self.simSubProcess.pid), signal.SIGTERM)
                except:
                    print("\n")
                self.simSubProcess = None
            self.mmaptx.cleanup()