# -*- coding: utf-8 -*-
# pNMS (ver4.0)

import pNMS.application_programming as cf
     
if __name__ == '__main__':

    ### [MODEL SELECTION]
    pool_count=10
    cf.createPool('musclefibers', pool_count)
    print ("Muscle-tendon fiber instances constructed")


    ### [MODEL PARAMETER SETTING]
    muscleFile='./parameters/MF_Parameters/MF_Parameters_3.0_SOL.csv'    
    cf.setParameters(muscleFile, init=True)
    print ("Model parameter values set")


    ### [SIMULATION CONDITION SETTING]
    cf.setSimulTimes(0, 2500, 0.05)
    print ("Simulation time set")
    
    # MTU state variables: CaSR, CaSRCS, CaSP, CaSPB, CaSPT, C1, C2, A, XCE
    muscle_ivalues = [0.0025, 0., 1.e-10, 0., 0., 0.154, 0.11, 0., -8]
    cf.setInitialValues(muscle_ivalues)
    print ("Initial value setting OK")
    
    
    ### [INPUT PARAMETER SETTING]
    # Fig S3 A
    cf.genSpikeSignals('Random', 100, 200, 1, 0) # twitch
    cf.setSpikeSignals(1, 1)  # (MF_#, ~ MF_#)
    cf.plotSpikeSignal()
    
    # # Fig S3 B
    cf.genSpikeSignals('Random', 100, 810, 10, 0) # 10hz
    cf.setSpikeSignals(2, 2)  # (MF_#, ~ MF_#)
    cf.plotSpikeSignal()
    
    # # Fig S3 C
    cf.genSpikeSignals('Random', 100, 810, 20, 0) # 20hz
    cf.setSpikeSignals(3, 3)  # (MF_#, ~ MF_#)
    cf.plotSpikeSignal()
    
    # # Fig S3 D
    cf.genSpikeSignals('Random', 100, 910, 40, 0) # 40hz
    cf.setSpikeSignals(4, 4)  # (MF_#, ~ MF_#)
    cf.plotSpikeSignal()
    
    # # Fig S3 A-D
    cf.genMuscleLengthSignals('Isometric', ivalue=-8)
    cf.setMuscleLengthSignals(1, 4)  # (MF_#, ~ MF_#)
    cf.plotMuscleLengthSignal()
    
    # Fig S3 E-H
    cf.genSpikeSignals('Random', 255, 806, 100, 0) # full excitation
    cf.setSpikeSignals(5, 8)  # (MF_#, ~ MF_#)
    cf.plotSpikeSignal()
    
    # Fig S3 E
    cf.genMuscleLengthSignals('Isokinetic', itime=616, ftime=730, ivalue=0, fvalue=-16) # shortennig
    cf.setMuscleLengthSignals(5, 5)  # (MF_#, ~ MF_#)
    cf.plotMuscleLengthSignal()
    
    # Fig S3 F
    cf.genMuscleLengthSignals('Isokinetic', itime=600, ftime=780, ivalue=-14, fvalue=0) # lengthning
    cf.setMuscleLengthSignals(6, 6)  # (MF_#, ~ MF_#)
    cf.plotMuscleLengthSignal()
    
    # # Fig S3 G
    cf.genMuscleLengthSignals('Isokinetic', ivalue=-8, fvalue=-8) # isometric
    cf.setMuscleLengthSignals(7, 7)  # (MF_#, ~ MF_#)
    cf.plotMuscleLengthSignal()
    
    # Fig S3 H
    cf.genMuscleLengthSignals('Isokinetic', itime=697, ftime=742, ivalue=-5, fvalue=-7) # brief step
    cf.setMuscleLengthSignals(8, 8)  # (MF_#, ~ MF_#)
    cf.plotMuscleLengthSignal()
    
    # Fig S3 I
    cf.genSpikeSignals('Random', 100, 1420, 10, 0) # 10hz
    cf.setSpikeSignals(9, 9)  # (MF_#, ~ MF_#)
    cf.plotSpikeSignal()
    
    # Fig S3 J
    cf.genSpikeSignals('Random', 100, 1510, 30, 0) # 30hz
    cf.setSpikeSignals(10, 10)  # (MF_#, ~ MF_#)
    cf.plotSpikeSignal()

    # Fig S3 I-J
    filePath = './parameters/Xm/locomotor-like movement Xm.csv'
    cf.importMuscleLengthSignals(filePath)
    cf.setMuscleLengthSignals(9, 10)  # (MF_#, ~ MF_#)
    cf.plotMuscleLengthSignal()
    
    
    ### [PARALLEL CONDITION SETTING]
    node_list = []
    cf.setComputeNode(node_list)
    print ("PP set")
    
    
    ### [PARALLEL SIMULATION]
    cf.runSimulation('autodetect', None)
    print ("Simulation done")
    
    
    ### [RESULT ONLINE DISPLAY]
    scope_list=['F', 'MUAP']
    cf.plotSimulResult(scope_list)
    print ("Online result plotting done")
    
    
    ### [RESULT SAVING]
    savePath='./results/'
    fileName='musclefibers'
    cf.saveSimulationResults(savePath, fileName)
    print ("Result saving done")
    
    
    ### [RESULT OFFLINE DISPLAY]
    dirPath = './results/'
    display = 'Combined'
    fileName = 'musclefibers'
    cf.plotImportData(dirPath, fileName, display, scope_list)
    print ("Offline result plotting done")