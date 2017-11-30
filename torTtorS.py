from random import randint
from random import gammavariate


def runSimulation(moleculesPerGen):
    
    #Initialize data structures (for 10 doublings)
    TorSLineage = [i for i in range(11)]
    TorTLineage = [i for i in range(11)]
    TorSTLineage = [i for i in range(11)]
    YFPLineage = [i for i in range(11)]
    
    for i in range(11):
        
        TorSLineage[i] = [0 for j in range(2**i)]
        TorTLineage[i] = [0 for j in range(2**i)]
        TorSTLineage[i] = [0 for j in range(2**i)]
        YFPLineage[i] = [0 for j in range(2**i)]    
    
    #Initialize the parent cell
        #Draw from gamma distribution to determine how much TorS and TorT are synthesized
        #Parameters for gamma distribution are alpha (bursts per generation) and beta (proteins produced per burst)
        #(See Friedman et al. 2006, Phys Rev Lett)
    TorSLineage[0][0] = int(round(gammavariate(moleculesPerGen,1)))
    TorTLineage[0][0] = int(round(gammavariate(moleculesPerGen,1)))
    TorRtot = 200.0 #Arbitrary
    
    #Make all TorS and TorT form TorST complex
    TorSTLineage[0][0] = min(TorSLineage[0][0], TorTLineage[0][0])
    TorSLineage[0][0] -= TorSTLineage[0][0]
    TorTLineage[0][0] -= TorSTLineage[0][0]
    
    #For parent cell assume YFP = TorR~P and equilibrium TorR~P levels
    if (TorSTLineage[0][0] != 0):
        YFPLineage[0][0] = TorSTLineage[0][0] * TorRtot / (TorSTLineage[0][0] + TorSLineage[0][0])
    else:
        YFPLineage[0][0] = 0
    
    
    #Loop through doublings
    for generation in range(10):
    
        #Loop through each cell in the population    
        for cell in range(2**generation):
    
            #Make all TorS and TorT form TorST complex (note: this is redundant for the parent cell [0][0])
            NewSTComplex = min(TorSLineage[generation][cell], TorTLineage[generation][cell])
            TorSLineage[generation][cell] -= NewSTComplex
            TorTLineage[generation][cell] -= NewSTComplex
            TorSTLineage[generation][cell] += NewSTComplex
    
            #TorR~P reaches its equilibrium value instantly
            if ((TorSTLineage[generation][cell] != 0) and (generation != 0) and (cell != 0)): #Parent cell is exempt
                TorRPhos = TorSTLineage[generation][cell] * TorRtot / (TorSTLineage[generation][cell] + TorSLineage[generation][cell])
            else:
                TorRPhos = 0        
    
            #Pull from gamma distribution to determine how much TorS and TorT are synthesized
            TorSProduced = int(round(gammavariate(moleculesPerGen,1)))
            TorTProduced = int(round(gammavariate(moleculesPerGen,1)))
            
            #Make all TorS and TorT form TorST complex
            NewSTComplex = min(TorSProduced, TorTProduced)
            TorSAfterGrowth = TorSLineage[generation][cell] + TorSProduced - NewSTComplex
            TorTAfterGrowth = TorTLineage[generation][cell] + TorTProduced - NewSTComplex
            TorSTAfterGrowth = TorSTLineage[generation][cell] + NewSTComplex
            
            #Synthesize 1 unit of YFP per unit of TorR~P                                       
            YFPAfterGrowth = YFPLineage[generation][cell] + TorRPhos 
    
            #Randomly assign TorS, TorT, and TorST to each daughter cell
            for molecule in range(TorSAfterGrowth):
    
                if randint(0, 1) == 0:
                    TorSLineage[generation + 1][cell * 2] += 1
                    
                else:                
                    TorSLineage[generation + 1][(cell * 2) + 1] += 1
    
            for molecule in range(TorTAfterGrowth):
    
                if randint(0, 1) == 0:
                    TorTLineage[generation + 1][cell * 2] += 1
                    
                else:                
                    TorTLineage[generation + 1][(cell * 2) + 1] += 1
    
            for molecule in range(TorSTAfterGrowth):
    
                if randint(0, 1) == 0:
                    TorSTLineage[generation + 1][cell * 2] += 1
                    
                else:                
                    TorSTLineage[generation + 1][(cell * 2) + 1] += 1
    
            #Divide YFP evenly between daughters               
            YFPLineage[generation + 1][cell * 2] = YFPAfterGrowth / 2
            YFPLineage[generation + 1][(cell * 2) + 1] = YFPAfterGrowth / 2
    
    #Return final generation's YFP values
    return YFPLineage[10]


def main():
    
    moleculesPerGen = float(raw_input("Enter number of TorS and TorT bursts per generation: "))
    numberOfRuns = int(raw_input("Enter number of times to run simulation: "))
    
    YFPFile = open('YFP_results.txt', 'w')
    YFPList = [0] * numberOfRuns #Initialize lines in output file
    
    for run in range(numberOfRuns):
        YFPList[run] = runSimulation(moleculesPerGen)
        
        for index in range(len(YFPList[run])):
            YFPFile.write(str(YFPList[run][index]) + '\n') #Write each final YFP value from all simulation runs in a column
        
    YFPFile.close()


main()