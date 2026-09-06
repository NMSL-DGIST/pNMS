# -*- coding: utf-8 -*-
# pNMS (ver4.0)

import pNMS.application_programming as cf

if __name__ == '__main__':
       

    ### [MODEL SELECTION]
    pool_count=10 #pool_count 10으로 변경
    cf.createPool('motoneuron', pool_count)
    print ("Motoneuron instances constructed")
    

    ### [MODEL PARAMETER SETTING]
    motorFile='./parameters/MN_Parameters/MN_Parameters_3.0_FR.csv'
    cf.setParameters(motorFile, init=True)
    print ("Model parameter values set")


    ### [SIMULATION CONDITION SETTING]
    cf.setSimulTimes(0, 20000, 0.1)
    print ("Simulation time set")
    
    # MN state variables: vs, sca, snam, snah, snapm, skdr, scam, scah, shm, vd,  dca, dcal, dnam, dnah, dnapm, dkdr, dkcam, dcam, dcah, dhm
    moto_ivalues=[-70., 0.0001, 0.001, 0.5829, 0.001, 0.1239, 0.004199, 0.9219, 0., -70., 0.0001, 0.001, 0.001, 0.5829, 0.001, 0.1239, 0.001, 0.004199, 0.9219, 0., ]
    cf.setInitialValues(moto_ivalues)
    print ("Initial values of state variables set")

    
    ### [INPUT PARAMETER SETTING]
     # Fig S1 A    
    # heav = [6.795, 0.5, 1000., 1300., -0.3, 3000., 3300., 0.37, 5000., 9000., 0.45, 6500., 6800., -0.37, 9000., 9300., 17.]
    # cf.genNeuronInputSignals('Step', heav_param=heav)
    # cf.setNeuronInputSignals(1, pool_count)
    # cf.plotNeuronInputSignal()
    # print ("Step Isoma set")

     # Fig S1 B
    # for i in range(1, pool_count+1):
    #    cf.genNeuronInputSignals('Ramp', 5000., 20., 0.)
    #    cf.setNeuronInputSignals(i, i)
    #    cf.plotNeuronInputSignal()
    # print ("Ramp Isoma set")

     # Fig S2 A
    # syn_heav_param = [0.029, 0.02, 1000., 1300., -0.028, 3000., 3300., 0.013, 5000., 9000., 0.03, 6500., 6800., -0.1,
    #                  9000., 9300., 1.]
    # for i in range(1, pool_count+1):
    #     cf.genSynConSignals('Dendrite', 'Excitatory', 'Step', syn_heav_param)
    #     cf.setSynConSignals('Dendrite', i, i)
    #     cf.plotSynConSignal('Dendrite')
    # print ("Step Iesyn set")

    # Fig S2 B
    syn_heav_param = []
    for i in range(1, pool_count+1):
        cf.genSynConSignals('Dendrite', 'Excitatory', 'Ramp', syn_heav_param, 0, 0.12, 10000, 0.5, 0.03, False)
        cf.setSynConSignals('Dendrite', i, i)
        cf.plotSynConSignal('Dendrite')
    print ("Ramp Iesyn set")

    
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
    

    ### [OFFLINE RESULT DISPLAY]
    dirPath = './results/'
    display = 'Combined'
    fileName = 'Motoneuron'
    cf.plotImportData(dirPath, fileName, display, scope_list)
    print ("Offline result display done")