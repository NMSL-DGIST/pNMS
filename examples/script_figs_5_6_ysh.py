# -*- coding: utf-8 -*-
# pNMS (ver4.0)

import pNMS.application_programming as cf

if __name__ == '__main__':
       
    ### [MODEL SELECTION]
    pool_count=10
    cf.createPool('motoneuron', pool_count)
    print ("Motoneuron instances constructed") ##print( ) 함수 괄호 추가  
    
    
    ### [MODEL PARAMETER SETTING]
    motorFile='./parameters/MN_Parameters/MN_Parameters_3.0_FR.csv'

    # default range model parameter (RMP) distribution
    RN_file='./parameters/RN Dpath/RN_decrease_exp.csv' # RN: 0.4 ~ 4.0 (MOhm)
    Dpath_func = 'linear' # f(RN): 0.1 ~ 1.1 (mm)
    Dpath_params = [0.25, 0.1] # slope, y-intercept
    
    cf.setParameters(motorFile, RN_file=RN_file, Dpath_func=Dpath_func, Dpath_params=Dpath_params)
    print ("Model parameter values set")
    
    ## Fig 5
 
    cf.plot_params('RN') 
    cf.plot_params('Dpath')
    cf.plot_params('gms')
    cf.plot_params('gmd')
    cf.plot_params('gc')
    cf.plot_params('cms')
    cf.plot_params('cmd')
    cf.plot_params('sf')
    cf.plot_params('sgna')
    cf.plot_params('dgcal')
    cf.plot_params('SNM')
    cf.plot_params('dgkca')    
    ### [SIMULATION CONDITION SETTING]
    cf.setSimulTimes(0, 20000, 0.05)
    print ("Simulation time set")

    moto_ivalues=[-70., 0.0001, 0.001, 0.5829, 0.001, 0.1239, 0.004199, 0.9219, 0., -70., 0.0001, 0.001, 0.001, 0.5829, 0.001, 0.1239, 0.001, 0.004199, 0.9219, 0., ]
    cf.setInitialValues(moto_ivalues)
    print ("Initial values of state variables set")
    
    ### [INPUT PARAMETER SETTING]
    
    # Fig 6 A: ramp current inection at the soma
    # for i in range(1, pool_count+1):
    #      cf.genNeuronInputSignals('Ramp', 5000., 20., 0.)
    #      cf.setNeuronInputSignals(i, i)
    #      cf.plotNeuronInputSignal()
    # print ("Ramp Isoma set")

    # Fig 6 B: ramp synaptic input at the dendrite
    syn_heav_param = []
    for i in range(1, pool_count+1):
        cf.genSynConSignals('Dendrite', 'Excitatory', 'Ramp', syn_heav_param, 0, 0.12, 10000, 0.5, 0.03, False)
        cf.genSynConSignals('Dendrite', 'Inhibitory', 'Ramp', syn_heav_param, 0, 0, 10000, 20, 0.035, False)
        cf.setSynConSignals('Dendrite', i, i)
        cf.plotSynConSignal('Dendrite')
    print ("Intracellular current inputs set")  
    

    

    
    ### [PARALLEL CONDITION SETTING]
    node_list=[]
    cf.setComputeNode(node_list)
    print ("Parallel simulation environment set")

    
    ### [PARALLEL SIMULATION]
    cf.runSimulation('autodetect',None)
    print ("Simulation done")
        
    
    ### [ONLINE RESULT DISPLAY]
    scope_list=['V_soma', 'Firing_rate']
    cf.plotSimulResult(scope_list)
    print ("Online result display done")


    ### [RESULT SAVING]
    savePath='./results/'
    fileName='Motoneuron'
    cf.saveSimulationResults(savePath, fileName)
    print ("Simulation results saved")
    

    ### [RESULT OFFLINE DISPLAY]
    dirPath = './results/'
    display = 'Combined'
    fileName = 'Motoneuron'
    cf.plotImportData(dirPath, fileName, display, scope_list)
    print ("Offline result display done")
