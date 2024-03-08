import subprocess
import mmapTx
import os
import signal

import sys,os


class SimCellPack:
# This class calls a cell simulation executable with the specified parameters
    # Must be called from within rl_balance
    
    def __init__(self, cellModel, numCells, simCycles, seed, balancing, sampleFactor, utilization):
        # Initialize the Cell Simulation class
        
        self.simSubProcess = None
        self.mmaptx = None

        self.numCells = numCells
        self.simCycles = simCycles

        self.seed = seed

        self.cellModels = ["A12", "ATL", "E1", "E2", "P14", "SAM"]
        self.balancingTypes = ["active", "passive"]

        try:
            index = os.getcwd().split("\\").index("rl_balance")
        except:
            exit("rl_balance folder not found, make sure to run this script inside rl_balance folder.")

        self.modelsDir = "\\".join(os.getcwd().split("\\")[:index+1]) + "\\cell_simulation"
        
        self.simExecutable = "\"" + self.modelsDir + "\\runPackSim.exe\""

        if cellModel not in self.cellModels:
            exit("Cell model specified not in the list: " + self.cellModels)
        else:
            self.cellModelPath = "\"" + self.modelsDir + "\\data\\cell_models\\" + cellModel + "model.mat\""

        if balancing not in self.balancingTypes:
            exit("Balancing type specified not in the list: " + self.balancingTypes)
        else:
            self.balancing = balancing
    
        self.profilePath = "\"" + self.modelsDir + "\\data\\drive_cycle_profiles\\us60Power.mat\""

        self.cellRandOpts = "[0,1,1,1,1,1]"

        self.sampleFactor = sampleFactor

        self.utilization = utilization
        
        
    def startSim(self):
    # Starts the execution of the cell executable simulator based on the configured parameters
        self.simCmd = 0
        
        self.cmdExe = " ".join([self.simExecutable, str(self.numCells), str(self.simCycles),        \
                                self.profilePath, self.cellModelPath, str(self.seed), self.balancing,    \
                                self.cellRandOpts, str(self.sampleFactor), str(self.utilization)])
                   
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