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

    # default range model parameter (RMP) distribution    
    p0_func = 'inc_convex' 
    p0_params = [0.01, 5]   # x range for target curvature
    p0_range = [0.1, 10]    # P0.5: 0.1 ~ 10 (N)
    
    # user-defined range model parameter (RMP) distribution
    tau1_func = 'linear'            # f(P0): 3.0 ~ 1.0
    tau1_params = [-0.2, 3.02]      # slope, y-intercept
    tau2_func = 'linear'            # f(P0): 25 ~ 13
    tau2_params = [-1.21, 25.12] 
    KSE_func = 'linear'             # f(P0): 0.4 ~ 0.16
    KSE_params = [-0.024, 0.402]
    AM_func = 'linear'              # f(P0): 0.1 ~ 0.5
    AM_params = [0.04, 0.096]
    LM_func = 'linear'              # f(P0): 1.0 ~ 0.5
    LM_params = [-0.051, 1.005]
    cv_func = 'linear'              # f(P0): 10 ~ 100
    cv_params = [9.09, 9.09]
    phi1_func = 'linear'            # f(P0): 0.03 ~ 0.0025
    phi1_params = [-0.003, 0.03]
    phi3_func = 'linear'            # f(P0): 0.01 ~ 0.000125
    phi3_params = [-0.001, 0.01]
    C1i_func = 'linear'             # f(P0): 0.12 ~ 0.154
    C1i_params = [0.004, 0.12]
    C1n1_func = 'linear'            # f(P0): 0 ~ 0.01
    C1n1_params = [0.001, 0]
    C1n4_func = 'linear'            # f(P0): 0.0001 ~ 85
    C1n4_params = [8.5858, -0.8585]
    C2i_func = 'linear'             # f(P0): 0.093 ~ 0.11
    C2i_params = [0.006, 0.089]
    C2n1_func = 'linear'            # f(P0): 0 ~ -0.0315
    C2n1_params = [-0.004, 0.0004]
    C2n4_func = 'linear'            # f(P0): 0.0001 ~ 70
    C2n4_params = [7.0707, -0.7070]
    alpha_i_func = 'linear'         # f(P0): 2 ~ 1.65
    alpha_i_params = [-0.04, 2.0]
    beta_func = 'linear'            # f(P0): 0.47 ~ 0.09
    beta_params = [-0.04, 0.5]
    gamma_func = 'linear'           # f(P0): 0.001 ~ 0
    gamma_params = [-0.0001, 0.001]
    g1_func = 'linear'              # f(P0): -8 ~ -0.8332
    g1_params = [0.73, -8.07]
    g2_func = 'linear'              # f(P0): 22 ~ 17.2907
    g2_params = [-0.51, 22.05]
    a0_func = 'linear'              # f(P0): 0.1 ~ 0.004
    a0_params = [-0.0097, 0.1010]
    b0_func = 'linear'              # f(P0): 24.35 ~ 99.7
    b0_params = [7.68, 23.23]
    c0_func = 'linear'              # f(P0): -0.32 ~ -0.58
    c0_params = [-0.0263, -0.3174]
    d0_func = 'linear'              # f(P0): 30.3 ~ 42.2
    d0_params = [1.31, 29.87]
    
    cf.setParameters(muscleFile, p0_func=p0_func, p0_params=p0_params, p0_range=p0_range,
                     KSE_func=KSE_func, KSE_params=KSE_params, 
                     cv_func=cv_func, cv_params=cv_params,
                     AM_func=AM_func, AM_params=AM_params, 
                     LM_func=LM_func, LM_params=LM_params, 
                     tau1_func=tau1_func, tau1_params=tau1_params, 
                     tau2_func=tau2_func, tau2_params=tau2_params,
                     phi1_func=phi1_func, phi1_params=phi1_params, 
                     phi3_func=phi3_func, phi3_params=phi3_params,
                     C1i_func=C1i_func, C1i_params=C1i_params, 
                     C1n1_func=C1n1_func, C1n1_params=C1n1_params, 
                     C1n4_func=C1n4_func, C1n4_params=C1n4_params,
                     C2i_func=C2i_func, C2i_params=C2i_params,
                     C2n1_func=C2n1_func, C2n1_params=C2n1_params,
                     C2n4_func=C2n4_func, C2n4_params=C2n4_params,
                     alpha_i_func=alpha_i_func, alpha_i_params=alpha_i_params, 
                     beta_func=beta_func, beta_params=beta_params, 
                     gamma_func=gamma_func, gamma_params=gamma_params,
                     g1_func=g1_func, g1_params=g1_params, 
                     g2_func=g2_func, g2_params=g2_params, 
                     a0_func=a0_func, a0_params=a0_params, 
                     b0_func=b0_func, b0_params=b0_params, 
                     c0_func=c0_func, c0_params=c0_params, 
                     d0_func=d0_func, d0_params=d0_params)
    print( "Model parameter values set")
        
    # Fig 7
    # cf.plot_params('p0')
    # cf.plot_params('KSE')
    # cf.plot_params('cv')
    # cf.plot_params('AM')
    # cf.plot_params('LM')
    # cf.plot_params('tau1')
    # cf.plot_params('tau2')
    # cf.plot_params('phi1')
    # cf.plot_params('phi3')
    # cf.plot_params('C1i')
    # cf.plot_params('C1n1')
    # cf.plot_params('C1n4')
    # cf.plot_params('C2i')
    # cf.plot_params('C2n1')
    # cf.plot_params('C2n4')
    # cf.plot_params('alpha_i')
    # cf.plot_params('beta')
    # cf.plot_params('gamma')
    # cf.plot_params('g1')
    # cf.plot_params('g2')
    # cf.plot_params('a0')
    # cf.plot_params('b0')
    # cf.plot_params('c0')
    # cf.plot_params('d0')
        
    
    ### [SIMULATION CONDITION SETTING]
    cf.setSimulTimes(0, 2000, 0.05)
    print ("Simulation time set")
    
    muscle_ivalues = [0.0025, 0., 1.e-10, 0., 0., 0.154, 0.11, 0., -8]
    cf.setInitialValues(muscle_ivalues)
    print ("Initial values set")


    ### [INPUT PARAMETER SETTING]
     #Fig 8 A
    cf.genSpikeSignals('Random', 100, 200, 1, 0) # twitch
    cf.setSpikeSignals(1, pool_count)  # (MF_#, ~ MF_#)
    cf.plotSpikeSignal()

    cf.genMuscleLengthSignals('Isometric', ivalue=-8)
    cf.setMuscleLengthSignals(1, pool_count)  # (MF_#, ~ MF_#)
    cf.plotMuscleLengthSignal()

     # Fig 8 B
    # cf.genSpikeSignals('Random', 100, 810, 10, 0) # 10hz
    # cf.setSpikeSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotSpikeSignal()

    # cf.genMuscleLengthSignals('Isometric', ivalue=-8)
    # cf.setMuscleLengthSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotMuscleLengthSignal()

     # Fig 8 C
    # cf.genSpikeSignals('Random', 100, 810, 20, 0) # 20hz
    # cf.setSpikeSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotSpikeSignal()

    # cf.genMuscleLengthSignals('Isometric', ivalue=-8)
    # cf.setMuscleLengthSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotMuscleLengthSignal()

     # Fig 8 D
    # cf.genSpikeSignals('Random', 100, 910, 40, 0) # 40hz
    # cf.setSpikeSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotSpikeSignal()

    # cf.genMuscleLengthSignals('Isometric', ivalue=-8)
    # cf.setMuscleLengthSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotMuscleLengthSignal()

     # Fig 8 E
    # cf.genSpikeSignals('Random', 255, 806, 100, 0) # full excitation
    # cf.setSpikeSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotSpikeSignal()

    # cf.genMuscleLengthSignals('Isokinetic', itime=616, ftime=730, ivalue=0, fvalue=-16) # shortennig
    # cf.setMuscleLengthSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotMuscleLengthSignal()

     # Fig 8 F
    # cf.genSpikeSignals('Random', 255, 806, 100, 0) # full excitation
    # cf.setSpikeSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotSpikeSignal()

    # cf.genMuscleLengthSignals('Isokinetic', itime=600, ftime=780, ivalue=-14, fvalue=0) # lengthning
    # cf.setMuscleLengthSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotMuscleLengthSignal()

     # Fig 8 G
    # cf.genSpikeSignals('Random', 255, 806, 100, 0) # full excitation
    # cf.setSpikeSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotSpikeSignal()
    
    # cf.genMuscleLengthSignals('Isokinetic', ivalue=-8, fvalue=-8) # isometric
    # cf.setMuscleLengthSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotMuscleLengthSignal()

     # Fig 8 H
    # cf.genSpikeSignals('Random', 255, 806, 100, 0) # full excitation
    # cf.setSpikeSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotSpikeSignal()
    
    # cf.genMuscleLengthSignals('Isokinetic', itime=697, ftime=742, ivalue=-5, fvalue=-7) # brief step
    # cf.setMuscleLengthSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotMuscleLengthSignal()

     # Fig 8 I
    # cf.genSpikeSignals('Random', 100, 1420, 10, 0) # 10hz
    # cf.setSpikeSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotSpikeSignal()
    
    # filePath = './parameters/Xm/locomotor-like movement Xm.csv'
    # cf.importMuscleLengthSignals(filePath)
    # cf.setMuscleLengthSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotMuscleLengthSignal()
    
    # Fig 8 J
    # cf.genSpikeSignals('Random', 100, 1510, 30, 0) # 30hz
    # cf.setSpikeSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotSpikeSignal()
    
    # filePath = './parameters/Xm/locomotor-like movement Xm.csv'
    # cf.importMuscleLengthSignals(filePath)
    # cf.setMuscleLengthSignals(1, pool_count)  # (MF_#, ~ MF_#)
    # cf.plotMuscleLengthSignal()
     
  
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
    cf.plotImportData(dirPath, fileName, 'Sum', scope_list)
    print( "Offline result plotting done")