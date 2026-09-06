# -*- coding: utf-8 -*-
# pNMS (ver4.0)

import pNMS.application_programming as cf

if __name__ == '__main__':
  
    ### [MODEL SELECTION]
    pool_count=10
    cf.createPool('motorunit', pool_count)
    print ("Motor unit instances constructed")
    
    
    ### [MODEL PARAMETER SETTING]
    motorFile='./parameters/MN_Parameters/MN_Parameters_3.0_FR.csv'
    muscleFile='./parameters/MF_Parameters/MF_Parameters_3.0_SOL.csv'
    
    
    # default range model parameter (RMP) distribution
    RN_file='./parameters/RN Dpath/RN_decrease_exp.csv' # RN: 0.4 ~ 4.0 (MOhm)
    Dpath_func = 'linear' # f(RN): 0.1 ~ 1.1 (mm)
    Dpath_params = [0, 0.6] # slope, y-intercept
    p0_func = 'linear' # f(RN): 10 ~ 0.1 (N)
    p0_params = [-2.75, 11.1] # slope, y-intercept
    
    # user-defined range model parameter (RMP) distribution
    dgkca_func = 'linear' # f(RN): 0.150 ~ 0.001 (mS/cm2) 
    dgkca_params = [0, 0.15] # slope, y-intercept

    cf.setParameters(motorFile, muscleFile, RN_file=RN_file, Dpath_func=Dpath_func, Dpath_params=Dpath_params, p0_func=p0_func, p0_params=p0_params, dgkca_func=dgkca_func, dgkca_params=dgkca_params)
    print ("Model parameter values set")
    
    # Fig 9
    # cf.plot_params('RN', unitType='motorunit')
    # cf.plot_params('Dpath')
    # cf.plot_params('gms')
    # cf.plot_params('gmd')
    # cf.plot_params('gc')
    # cf.plot_params('cms')
    # cf.plot_params('cmd')
    # cf.plot_params('sf')
    # cf.plot_params('sgna')
    # cf.plot_params('dgcal')
    # cf.plot_params('SNM')
    # cf.plot_params('dgkca')
    # cf.plot_params('p0', unitType='motorunit')
    # cf.plot_params('tau1')
    # cf.plot_params('tau2')
    # cf.plot_params('KSE')
    # cf.plot_params('AM')
    # cf.plot_params('LM')
    # cf.plot_params('cv')
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
    cf.setSimulTimes(0, 20000, 0.1)
    print ("Simulation time set")
    
    moto_ivalues=[-70., 0.0001, 0.001, 0.5829, 0.001, 0.1239, 0.004199, 0.9219, 0., -70., 0.0001, 0.001, 0.001, 0.5829, 0.001, 0.1239, 0.001, 0.004199, 0.9219, 0., ]
    muscle_ivalues = [0.0025, 0., 1.e-10, 0., 0., 0.154, 0.11, 0., -8]
    cf.setInitialValues(moto_ivalues, muscle_ivalues)  
    print ("Initial values set")


    ### [INPUT PARAMETER SETTING] 
    # Fig 10 A
    cf.genNeuronInputSignals('Ramp', 5000., 20., 0.)
    cf.setNeuronInputSignals(1, pool_count)
    cf.plotNeuronInputSignal()
    print ("Ramp Isoma set")
    
    cf.genMuscleLengthSignals('Isometric', ivalue=-8)
    cf.setMuscleLengthSignals(1, pool_count)
    cf.plotMuscleLengthSignal()
    print ("Muscle-tendon length set")

    # Fig 10 B
    # syn_heav_param = []
    # for i in range(1, pool_count+1):
    #    cf.genSynConSignals('Dendrite', 'Excitatory', 'Ramp', syn_heav_param, 0, 0.3, 10000, 20, 0.06, True)
    #    cf.genSynConSignals('Dendrite', 'Inhibitory', 'Ramp', syn_heav_param, 0.02, 0.02, 10000, 20, 0.035, True)
    #    cf.setSynConSignals('Dendrite', i, i)
    #    cf.plotSynConSignal('Dendrite')
    # print ("Ramp Iesyn and background Iisyn set")

    # cf.genMuscleLengthSignals('Isometric', ivalue=-8)
    # cf.setMuscleLengthSignals(1, pool_count)
    # cf.plotMuscleLengthSignal()
    # print ("Muscle length set")    
  
    # Fig 10 C
    # syn_heav_param =[]
    # for i in range(1, pool_count + 1):
    #     cf.genSynConSignals('Dendrite', 'Excitatory', 'Ramp', syn_heav_param, 0, 0.3, 10000, 20, 0.06, True)
    #     cf.genSynConSignals('Dendrite', 'Inhibitory', 'Ramp', syn_heav_param, 0.02, 0.165, 10000, 20, 0.035, True)
    #     cf.setSynConSignals('Dendrite', i, i) # motoneuron#1 to motoneuron#10
    #     cf.plotSynConSignal('Dendrite')
    # print( "Balanced Iesyn and Iisyn set")
    
    # cf.genMuscleLengthSignals('Isometric', ivalue=-8)
    # cf.setMuscleLengthSignals(1, pool_count)
    # cf.plotMuscleLengthSignal()
    # print ("Muscle length set")

    # Fig 10 D
    # syn_heav_param = []
    # for i in range(1, pool_count + 1):
    #     cf.genSynConSignals('Dendrite', 'Excitatory', 'Ramp', syn_heav_param, 0, 0.3, 10000, 20, 0.06, True)
    #     cf.genSynConSignals('Dendrite', 'Inhibitory', 'Ramp', syn_heav_param, 0.165, 0.02, 10000, 20, 0.035, True)
    #     cf.setSynConSignals('Dendrite', i, i) # motoneuron#16 to motoneuron#20
    #     cf.plotSynConSignal('Dendrite')
    # print ("Push-pull Iesyn and Iisyn set")
    
    # cf.genMuscleLengthSignals('Isometric', ivalue=-8)
    # cf.setMuscleLengthSignals(1, pool_count)
    # cf.plotMuscleLengthSignal()
    # print ("Muscle length set")
    
        
    ### [PARALLEL CONDITION SETTING]
    node_list=[]
    cf.setComputeNode(node_list)
    print ("PP set")

    
    ### [PARALLEL SIMULATION]
    cf.runSimulation('autodetect',None)
    print ("Simulation done")

    
    ### [RESULT ONLINE DISPLAY]
    scope_list=['V_soma', 'Firing_rate', 'F']
    cf.plotSimulResult(scope_list)
    print ("Online result plotting done")

    
    ### [RESULT SAVING]
    savePath='./results/'
    fileName='Motorunit'
    cf.saveSimulationResults(savePath, fileName)
    print ("Result saving done")
    
    
    ### [RESULT OFFLINE DISPLAY]
    dirPath = './results/'
    display = 'Combined'
    scope_list = ['V_soma', 'Firing_rate', 'F', 'MUAP']
    fileName = 'Motorunit'
    cf.plotImportData(dirPath, fileName, display, scope_list)
    
    display = 'Sum'
    scope_list = ['F', 'MUAP']
    fileName = 'Motorunit'
    cf.plotImportData(dirPath, fileName, display, scope_list)
    print ("Offline result plotting done")