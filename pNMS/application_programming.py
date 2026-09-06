# -*- coding: utf-8 -*-
# pMUPS (ver4.0)


import pNMS.simulation_engine as sim
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#v4.0
# import matplotlib as mpl # CHJ
# from cycler import cycler # CHJ
# mpl.rcParams['axes.prop_cycle'] = cycler('color', ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'purple', 'pink', 'brown']) # CHJ

import csv


## Model parameter files
motorFile='./parameters/MN_Parameters/MN_Parameters_3.0_FR.csv'
muscleFile='./parameters/MF_Parameters/MF_Parameters_3.0_SOL.csv'
# muscleFile='./parameters/MF_Parameters/MF_Parameters_3.0_MG.csv'

## Initial values for model state variables
# motoneuron: vs, sca, snam, snah, snapm, skdr, scam, scah, shm, vd,  dca, dcal, dnam, dnah, dnapm, dkdr, dkcam, dcam, dcah, dhm
MN_ivalue=[-70., 0.0001, 0.001, 0.5829, 0.001, 0.1239, 0.004199, 0.9219, 0., -70., 0.0001, 0.001, 0.001, 0.5829, 0.001, 0.1239, 0.001, 0.004199, 0.9219, 0., ]
# muscle-tendon unit: CaSR, CaSRCS, CaSP, CaSPB, CaSPT, C1, C2, A, XCE
MF_ivalue=[0.0025, 0., 1.e-7, 0., 0., 0.154, 0.11, 0., -8]

## Default coefficient values for heavian function
# coefficients: i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s
Is_heav_param=[6.795, 0.5, 1000., 1300., -0.3, 3000., 3300., 0.37, 5000., 9000., 0.45, 6500., 6800., -0.37, 9000., 9300., 17.]
Syn_heav_param=[0.029, 0.02, 1000.,1300., -0.028, 3000.,3300., 0.013, 5000.,9000., 0.03, 6500.,6800., -0.1, 9000.,9300., 1.]

## Homogeneous population modeling
def createPool(pool_type='motorunit', num=10):
    
    global MP, NSG, s_SCSG, d_SCSG, SSG, MSG
    
    if(pool_type=='motoneuron'):
        MP=sim.MotoneuronPool(unit_count=num)
        setParameters(motorFile, init=True)
        setInitialValues(MN_ivalue)
        setSimulTimes()
        NSG=sim.NeuronSignalGenerator(MP, num)
        s_SCSG = sim.SynConSignalGenerator(MP, num, 'Soma')
        d_SCSG = sim.SynConSignalGenerator(MP, num, 'Dendrite')

    elif(pool_type=='musclefibers'):
        MP=sim.MusclefibersPool(unit_count=num)
        setParameters(muscleFile, init=True)
        setInitialValues(MF_ivalue)
        setSimulTimes()
        SSG=sim.SpikeSignalGenerator(MP, num)
        MSG=sim.MuscleLengthSignalGenerator(MP, num)

    elif(pool_type=='motorunit'):
        MP=sim.MotorunitPool(unit_count=num)
        setParameters(motorFile, muscleFile, init=True)
        setInitialValues(MN_ivalue, MF_ivalue)
        setSimulTimes()
        NSG=sim.NeuronSignalGenerator(MP, num)
        s_SCSG = sim.SynConSignalGenerator(MP, num, 'Soma')
        d_SCSG = sim.SynConSignalGenerator(MP, num, 'Dendrite')
        MSG=sim.MuscleLengthSignalGenerator(MP, num)


## Heterogeneous population modeling
# 241210 and 241211
#def setParameters(param_file, MF_file=None, RN_file=None, Dpath_file=None, RN_func=None, Dpath_func=None, RN_params=None, Dpath_params=None, 
#                  p0_file=None, tau1_file=None,tau2_file=None,KSE_file=None,AM_file=None,LM_file=None, cv_file=None,
#                  phi1_file=None, phi3_file=None, C1i_file=None, C1n1_file=None, C1n4_file=None, C2i_file=None, C2n1_file=None, C2n4_file=None,
#                  C3_file=None, C4_file=None, C5_file=None, alpha_i_file=None, beta_file=None, gamma_file=None, g1_file=None, g2_file=None, a0_file=None, b0_file=None, c0_file=None, d0_file=None,
#                  p0_func=None, tau1_func=None,tau2_func=None,KSE_func=None,AM_func=None,LM_func=None, cv_func=None,
#                  phi1_func=None, phi3_func=None, C1i_func=None, C1n1_func=None, C1n4_func=None, C2i_func=None, C2n1_func=None, C2n4_func=None,
#                  C3_func=None, C4_func=None, C5_func=None, alpha_i_func=None, beta_func=None, gamma_func=None, g1_func=None, g2_func=None, a0_func=None, b0_func=None, c0_func=None, d0_func=None,
#                  p0_params=None, tau1_params=None,tau2_params=None,KSE_params=None,AM_params=None,LM_params=None, cv_params=None,
#                  phi1_params=None, phi3_params=None, C1i_params=None, C1n1_params=None, C1n4_params=None, C2i_params=None, C2n1_params=None, C2n4_params=None,
#                  C3_params=None, C4_params=None, C5_params=None, alpha_i_params=None, beta_params=None, gamma_params=None, g1_params=None, g2_params=None, a0_params=None, b0_params=None, c0_params=None, d0_params=None,
#                  init=False):

def setParameters(param_file, MF_file=None,
                  RN_file=None, Dpath_file=None, gms_file=None, gmd_file=None, gc_file=None, cms_file=None, cmd_file=None, sf_file=None, sgna_file=None, dgcal_file=None, dgkca_file=None, SNM_file=None,
                  p0_file=None, tau1_file=None,tau2_file=None,KSE_file=None,AM_file=None,LM_file=None, cv_file=None,
                  phi1_file=None, phi3_file=None, C1i_file=None, C1n1_file=None, C1n4_file=None, C2i_file=None, C2n1_file=None, C2n4_file=None,
                  C3_file=None, C4_file=None, C5_file=None, alpha_i_file=None, beta_file=None, gamma_file=None, g1_file=None, g2_file=None, a0_file=None, b0_file=None, c0_file=None, d0_file=None,
                  RN_func=None, Dpath_func=None, gms_func=None, gmd_func=None, gc_func=None, cms_func=None, cmd_func=None, sf_func=None, sgna_func=None, dgcal_func=None, dgkca_func=None, SNM_func=None,
                  p0_func=None, tau1_func=None,tau2_func=None,KSE_func=None,AM_func=None,LM_func=None, cv_func=None,
                  phi1_func=None, phi3_func=None, C1i_func=None, C1n1_func=None, C1n4_func=None, C2i_func=None, C2n1_func=None, C2n4_func=None,
                  C3_func=None, C4_func=None, C5_func=None, alpha_i_func=None, beta_func=None, gamma_func=None, g1_func=None, g2_func=None, a0_func=None, b0_func=None, c0_func=None, d0_func=None,
                  RN_params=None, Dpath_params=None, gms_params=None, gmd_params=None, gc_params=None, cms_params=None, cmd_params=None, sf_params=None, sgna_params=None, dgcal_params=None, dgkca_params=None, SNM_params=None,
                  p0_params=None, tau1_params=None,tau2_params=None,KSE_params=None,AM_params=None,LM_params=None, cv_params=None,
                  phi1_params=None, phi3_params=None, C1i_params=None, C1n1_params=None, C1n4_params=None, C2i_params=None, C2n1_params=None, C2n4_params=None,
                  C3_params=None, C4_params=None, C5_params=None, alpha_i_params=None, beta_params=None, gamma_params=None, g1_params=None, g2_params=None, a0_params=None, b0_params=None, c0_params=None, d0_params=None,
                  RN_range=None, Dpath_range=None, gms_range=None, gmd_range=None, gc_range=None, cms_range=None, cmd_range=None, sf_range=None, sgna_range=None, dgcal_range=None, dgkca_range=None, SNM_range=None,
                  p0_range=None, tau1_range=None,tau2_range=None,KSE_range=None,AM_range=None,LM_range=None, cv_range=None,
                  phi1_range=None, phi3_range=None, C1i_range=None, C1n1_range=None, C1n4_range=None, C2i_range=None, C2n1_range=None, C2n4_range=None,
                  C3_range=None, C4_range=None, C5_range=None, alpha_i_range=None, beta_range=None, gamma_range=None, g1_range=None, g2_range=None, a0_range=None, b0_range=None, c0_range=None, d0_range=None,
                  init=False):
                                            
    global MP, PG
    
    PG=sim.ParametersGenerator()
    unitType=MP.getUnitType()
    unitCount=MP.getUnitCount()
    index = 'Parameter'
    col = 'Value'
    
    if(unitType=='motoneuron' or unitType=='motorunit'):
        MN_file = param_file
        with open(MN_file, 'r') as file:
            csv_reader = csv.reader(file)
            data = [{index:param, col:value} for param, value in csv_reader]
        params_table = pd.DataFrame(data)[3:].set_index(index).astype('float') # first 3 rows skiped from Dataframe, str(csv_reader) -> float
        
        params = params_table['svl':][col].values
        mn_params=[]
        # tm = params_table.get_value('tm', col) # CHJ
        tm = params_table.loc['tm', col]
        # p = params_table.get_value('parea', col) # CHJ
        p = params_table.loc['parea', col]
        # dgcal scaling with neuromodulation factor
        # dgcal = params_table.get_value('SNM', col) * params_table.get_value('dgcal', col) # CHJ
        dgcal = params_table.loc['SNM', col] * params_table.loc['dgcal', col]
        # params_table.set_value('dgcal', col, dgcal) # CHJ
        params_table.loc['dgcal', col] = dgcal

        # homogeneous motoneuron population
        if(init==True):
            # RN = params_table.get_value('rn', col) # CHJ
            RN = params_table.loc['rn', col]
            # va_sd_dc = params_table.get_value('VAsdDC', col) # CHJ
            va_sd_dc = params_table.loc['VAsdDC', col]
            # va_sd_ac = params_table.get_value('VAsdAC', col) # CHJ
            va_sd_ac = params_table.loc['VAsdAC', col]
            # va_ds_dc = params_table.get_value('VAdsDC', col) # CHJ
            va_ds_dc = params_table.loc['VAdsDC', col]
            
            VA = (va_sd_dc, va_sd_ac, va_ds_dc)
            gms,gmd,gc,cmd,cms = PG.calCableParameters(RN, p, tm=tm, VA=VA)
            cable_params = np.array([gms,gmd,gc,cmd,cms])
            params = np.append(cable_params, params) # params = cable_params + params
            for i in range(unitCount):                
                mn_params.append(params)  # ndarray
        
        # heterogeneous motoneuron population
        else:
            if(RN_file != None): # RN from file
                RN_arr = PG.importParameters('RN', RN_file)
                table_x = np.arange(1, len(RN_arr) + 1)
                table_y = RN_arr
                interp_x = np.linspace(1, len(RN_arr), unitCount)
                interp_y = PG.cal_interp(table_x, table_y, 1, interp_x)
                RN_arr = interp_y
                PG.RN = RN_arr
                # print 'RN: ' # CHJ
                print('RN: ')
                # print RN_arr # CHJ
                print(RN_arr)
            else: # RN from function
                MN_arr = np.arange(1, unitCount+1)
                # 241210                
                #RN_arr = PG.genParameters('RN', RN_func, MN_arr, RN_params, unitCount)
                RN_arr = PG.genParameters('RN', RN_func, MN_arr, RN_params, unitCount, RN_range)

            if(Dpath_file != None): # Dpath from file
                Dpath_arr = PG.importParameters('Dpath', Dpath_file)
            else: # Dpath from function
                # 241210                
                #Dpath_arr = PG.genParameters('Dpath', Dpath_func, RN_arr, Dpath_params, unitCount)
                Dpath_arr = PG.genParameters('Dpath', Dpath_func, RN_arr, Dpath_params, unitCount, Dpath_range)
            # print 'Dpath: ' # CHJ
            print('Dpath: ')
            # print Dpath_arr # CHJ
            print(Dpath_arr)
            Dpath_arr[Dpath_arr < 0.02] = 0.02 # minimum Dpath (mm)
            # print Dpath_arr # CHJ
            print(Dpath_arr)
            
            # default valeus of range parameters
            gms_arr,gmd_arr,gc_arr,cmd_arr,cms_arr = PG.calCableParameters(RN_arr, p, Dpath=Dpath_arr) # inverse equations: f(RN, Dpath)
            sf_arr,sgna_arr,dgcal_arr,dgkca_arr = PG.interpParameters('motoneuron', RN_arr, Dpath_arr) # lookup tables: f(RN, Dpath)
            
            # 241211            
            # user-defined values of range parameters
            if (gms_file != None): # gms from file
                gms_arr = PG.importParameters('gms', gms_file)
            elif (gms_func != None): # gms from function
                gms_arr = PG.genParameters('gms', gms_func, RN_arr, gms_params, unitCount, gms_range) # f(RN)
            
            if (gmd_file != None): # gmd from file
                gmd_arr = PG.importParameters('gmd', gmd_file)
            elif (gmd_func != None): # gmd from function
                gmd_arr = PG.genParameters('gmd', gmd_func, RN_arr, gmd_params, unitCount, gmd_range) # f(RN)
            
            if (gc_file != None): # gc from file
                gc_arr = PG.importParameters('gc', gc_file)
            elif (gc_func != None): # gc from function
                gc_arr = PG.genParameters('gc', gc_func, RN_arr, gc_params, unitCount, gc_range) # f(RN)
            
            if (cms_file != None): # cms from file
                cms_arr = PG.importParameters('cms', cms_file)
            elif (cms_func != None): # cms from function
                cms_arr = PG.genParameters('cms', cms_func, RN_arr, cms_params, unitCount, cms_range) # f(RN)
            
            if (cmd_file != None): # cmd from file
                cmd_arr = PG.importParameters('cmd', cmd_file)
            elif (cmd_func != None): # cmd from function
                cmd_arr = PG.genParameters('cmd', cmd_func, RN_arr, cmd_params, unitCount, cmd_range) # f(RN)
                        
            if (sf_file != None): # sf from file
                sf_arr = PG.importParameters('sf', sf_file)
            elif (sf_func != None): # sf from function
                sf_arr = PG.genParameters('sf', sf_func, RN_arr, sf_params, unitCount, sf_range) # f(RN)
            
            if (sgna_file != None): # sgna from file
                sgna_arr = PG.importParameters('sgna', sgna_file)
            elif (sgna_func != None): # sgna from function
                sgna_arr = PG.genParameters('sgna', sgna_func, RN_arr, sgna_params, unitCount, sgna_range) # f(RN)
            
            if (dgcal_file != None): # dgcal from file
                dgcal_arr = PG.importParameters('dgcal', dgcal_file)
            elif (dgcal_func != None): # dgcal from function
                dgcal_arr = PG.genParameters('dgcal', dgcal_func, RN_arr, dgcal_params, unitCount, dgcal_range) # f(RN)
            
            if (dgkca_file != None): # dgkca from file
                dgkca_arr = PG.importParameters('dgkca', dgkca_file)
            elif (dgkca_func != None): # dgkca from function
                dgkca_arr = PG.genParameters('dgkca', dgkca_func, RN_arr, dgkca_params, unitCount, dgkca_range) # f(RN)
            
            if (SNM_file != None): # SNM from file
                SNM_arr = PG.importParameters('SNM', SNM_file)
                PG.dgcal = SNM_arr * dgcal_arr
            elif (SNM_func != None): # SNM from function
                SNM_arr = PG.genParameters('SNM', SNM_func, RN_arr, SNM_params, unitCount, SNM_range) # f(RN)
                PG.dgcal = SNM_arr * dgcal_arr
            # 241211
            
            num = 1
            for gms,gmd,gc,cmd,cms,sf,sgna,dgcal,dgkca in zip(gms_arr,gmd_arr,gc_arr,cmd_arr,cms_arr,sf_arr,sgna_arr,dgcal_arr,dgkca_arr):
                # params_table.set_value('sf', col, sf) # CHJ
                params_table.at['sf', col] = sf
                # params_table.set_value('df', col, sf) # df = sf # CHJ
                params_table.at['df', col] = sf
                # params_table.set_value('sgna', col, sgna) # CHJ
                params_table.at['sgna', col] = sgna
                # params_table.set_value('dgcal', col, dgcal) # CHJ
                params_table.at['dgcal', col] = dgcal
                # params_table.set_value('dgkca', col, dgkca) # CHJ
                params_table.at['dgkca', col] = dgkca
                cable_params = np.array([gms,gmd,gc,cmd,cms])
                whole_params = np.append(cable_params, params)  # ndarray: cable_params + params 
                mn_params.append(whole_params)  # ndarray
                
                # display of range parameter values
                # print 'MN %d : gms = %.3f, gmd = %.3f, gc = %.3f, cms = %.3f, cmd = %.3f, sgna = %.3f, sf = %.3f, dgcal = %.3f, dgkca = %.3f, df = %.3f' 
                # % (num,gms,gmd,gc,cms,cmd,sgna,sf,dgcal,dgkca, params_table.get_value('df', col)) # CHJ
                print('MN %d : gms = %.3f, gmd = %.3f, gc = %.3f, cms = %.3f, cmd = %.3f, sgna = %.3f, sf = %.3f, dgcal = %.3f, dgkca = %.3f, df = %.3f' 
                % (num, gms, gmd, gc, cms, cmd, sgna, sf, dgcal, dgkca, params_table.at['df', col]))
                num += 1

                
        if(unitType=='motoneuron'):
            MP.setParameters(mn_params)
    
    if(unitType=='musclefibers' or unitType=='motorunit'):
        if(unitType=='musclefibers'):
            MF_file = param_file
        with open(MF_file, 'r') as file:
            csv_reader = csv.reader(file)
            data = [{index:param, col:value} for param, value in csv_reader] # list of dictionary
        params_table = pd.DataFrame(data)[1:].set_index(index).astype('float') # first row skiped from Dataframe, str(csv_reader) -> float    
        params = params_table[:][col].values
        mf_params=[]    
        
        # homogeneous muscle-tendon unit population or motor unit population
        if(init==True):
            for i in range(unitCount):
                mf_params.append(params)

        # heterogeneous muscle-tendon unit population or motor unit population
#        else:
#            if(unitType=='musclefibers'):
#                MF_arr = np.arange(1, unitCount+1)
#
#                if(p0_file != None): # p0 from file
#                    p0_arr = PG.importParameters('p0', p0_file)
#                else: # p0 from function
#                    p0_arr = PG.genParameters('p0', p0_func, MF_arr, p0_params, unitCount)    # f(MF#)
#                
#                if(tau1_file != None): # tau1 from file
#                    tau1_arr = PG.importParameters('tau1', tau1_file)
#                else: # tau1 from function
#                    tau1_arr = PG.genParameters('tau1', tau1_func, p0_arr, tau1_params, unitCount) # f(p0)
#                
#                if(tau2_file != None): # tau2 from file
#                    tau2_arr = PG.importParameters('tau2', tau2_file)
#                else: # tau2 from function
#                    tau2_arr = PG.genParameters('tau2', tau2_func, p0_arr, tau2_params, unitCount) # f(p0)
#                
#                if(KSE_file != None): # KSE from file
#                    KSE_arr = PG.importParameters('KSE', KSE_file)
#                else: # KSE from function
#                    KSE_arr = PG.genParameters('KSE', KSE_func, p0_arr, KSE_params, unitCount) # f(p0)
#                
#                if(AM_file != None): # AM from file
#                    AM_arr = PG.importParameters('AM', AM_file)
#                else: # AM from function
#                    AM_arr = PG.genParameters('AM', AM_func, p0_arr, AM_params, unitCount) # f(p0)
#                
#                if(LM_file != None): # LM from file
#                    LM_arr = PG.importParameters('LM', LM_file)
#                else: # LM from function
#                    LM_arr = PG.genParameters('LM', LM_func, p0_arr, LM_params, unitCount) # f(p0)
#                
#                if(cv_file != None): # cv from file
#                    cv_arr = PG.importParameters('cv', cv_file)
#                else: # cv from function
#                    cv_arr = PG.genParameters('cv', cv_func, p0_arr, cv_params, unitCount) # f(p0)
#
#                if (phi1_file != None): # phi1 from file
#                    phi1_arr = PG.importParameters('phi1', phi1_file)
#                else: # phi1 from function
#                    phi1_arr = PG.genParameters('phi1', phi1_func, p0_arr, phi1_params, unitCount)  # f(p0)
#                
#                if (phi3_file != None): # phi3 from file
#                    phi3_arr = PG.importParameters('phi3', phi3_file)
#                else: # phi3 from function
#                    phi3_arr = PG.genParameters('phi3', phi3_func, p0_arr, phi3_params, unitCount)  # f(p0)
#
#                if (C1i_file != None): # C1i from file
#                    C1i_arr = PG.importParameters('C1i', C1i_file)
#                else: # C1i from function
#                    C1i_arr = PG.genParameters('C1i', C1i_func, p0_arr, C1i_params, unitCount)  # f(p0)
#                
#                if (C1n1_file != None): # C1n1 from file
#                    C1n1_arr = PG.importParameters('C1n1', C1n1_file)
#                else: # C1n1 from function
#                    C1n1_arr = PG.genParameters('C1n1', C1n1_func, p0_arr, C1n1_params, unitCount)  # f(p0)
#                
#                if (C1n4_file != None): # C1n4 from file
#                    C1n4_arr = PG.importParameters('C1n4', C1n4_file)
#                else: # C1n4 from function
#                    C1n4_arr = PG.genParameters('C1n4', C1n4_func, p0_arr, C1n4_params, unitCount)  # f(p0)
#
#                if (C2i_file != None): # C2i from file
#                    C2i_arr = PG.importParameters('C2i', C2i_file)
#                else: # C2i from function
#                    C2i_arr = PG.genParameters('C2i', C2i_func, p0_arr, C2i_params, unitCount)  # f(p0)
#                
#                if (C2n1_file != None): # C2n1 from file
#                    C2n1_arr = PG.importParameters('C2n1', C2n1_file)
#                else: # C2n1 from function
#                    C2n1_arr = PG.genParameters('C2n1', C2n1_func, p0_arr, C2n1_params, unitCount)  # f(p0)
#                
#                if (C2n4_file != None): # C2n4 from file
#                    C2n4_arr = PG.importParameters('C2n4', C1n4_file)
#                else: # C2n4 from function
#                    C2n4_arr = PG.genParameters('C2n4', C2n4_func, p0_arr, C2n4_params, unitCount)  # f(p0)
#
#                if (C3_file != None): # C3 from file
#                    C3_arr = PG.importParameters('C3', C3_file)
#                else: # C3 from function
#                    C3_arr = PG.genParameters('C3', C3_func, p0_arr, C3_params, unitCount)  # f(p0)
#                
#                if (C4_file != None): # C4 from file
#                    C4_arr = PG.importParameters('C4', C4_file)
#                else: # C4 from function
#                    C4_arr = PG.genParameters('C4', C4_func, p0_arr, C4_params, unitCount)  # f(p0)
#                
#                if (C5_file != None): # C5 from file
#                    C5_arr = PG.importParameters('C5', C5_file)
#                else: # C5 from function
#                    C5_arr = PG.genParameters('C5', C5_func, p0_arr, C5_params, unitCount)  # f(p0)
#
#                if (alpha_i_file != None): # alpha_i from file
#                    alpha_i_arr = PG.importParameters('alpha_i', alpha_i_file)
#                else: # alpha_i from function
#                    alpha_i_arr = PG.genParameters('alpha_i', alpha_i_func, p0_arr, alpha_i_params, unitCount)  # f(p0)
#                
#                if (beta_file != None): # beta from file
#                    beta_arr = PG.importParameters('beta', beta_file)
#                else: # beta from function
#                    beta_arr = PG.genParameters('beta', beta_func, p0_arr, beta_params, unitCount)  # f(p0)
#                
#                if (gamma_file != None): # gamma from file
#                    gamma_arr = PG.importParameters('gamma', gamma_file)
#                else: # gamma from function
#                    gamma_arr = PG.genParameters('gamma', gamma_func, p0_arr, gamma_params, unitCount)  # f(p0)
#
#                if (g1_file != None): # g1 from file
#                    g1_arr = PG.importParameters('g1', beta_file)
#                else: # g1 from function
#                    g1_arr = PG.genParameters('g1', g1_func, p0_arr, g1_params, unitCount)  # f(p0)
#                
#                if (g2_file != None): # g2 from file
#                    g2_arr = PG.importParameters('g2', g2_file)
#                else: # g2 from function
#                    g2_arr = PG.genParameters('g2', g2_func, p0_arr, g2_params, unitCount)  # f(p0)
#                
#                if (a0_file != None): # a0 from file
#                    a0_arr = PG.importParameters('a0', a0_file)
#                else: # a0 from function
#                    a0_arr = PG.genParameters('a0', a0_func, p0_arr, a0_params, unitCount)  # f(p0)
#
#                if (b0_file != None): # b0 from file
#                    b0_arr = PG.importParameters('b0', b0_file)
#                else: # b0 from function
#                    b0_arr = PG.genParameters('b0', b0_func, p0_arr, b0_params, unitCount)  # f(p0)
#
#                if (c0_file != None): # c0 from file
#                    c0_arr = PG.importParameters('c0', c0_file)
#                else: # c0 from function
#                    c0_arr = PG.genParameters('c0', c0_func, p0_arr, c0_params, unitCount)  # f(p0)
#
#                if (d0_file != None): # d0 from file
#                    d0_arr = PG.importParameters('d0', d0_file)
#                else: # d0 from function
#                    d0_arr = PG.genParameters('d0', d0_func, p0_arr, d0_params, unitCount)  # f(p0)
#            
#            elif(unitType=='motorunit'):
#                # 241208                
#                MF_arr = np.arange(1, unitCount+1)
#                if(p0_file != None): # p0 from file
#                    p0_arr = PG.importParameters('p0', p0_file)
#                else: # p0 from function
#                    p0_arr = PG.genParameters('p0', p0_func, RN_arr, p0_params, unitCount)    # f(MF#)
#                # print 'p0: ' # CHJ
#                print('p0: ')
#                # print p0_arr # CHJ
#                print(p0_arr)
#                
#                # 241208
#                # interpolation using lookup tables with linear relationship
#                #p0_arr,tau1_arr,tau2_arr,KSE_arr,AM_arr,LM_arr,cv_arr,phi1_arr,phi3_arr,\
#                #C1i_arr,C1n1_arr,C1n4_arr,C2i_arr,C2n1_arr,C2n4_arr,C3_arr,C4_arr,C5_arr,\
#                #alpha_i_arr,beta_arr,gamma_arr,g1_arr,g2_arr,a0_arr,b0_arr,c0_arr,d0_arr = PG.interpParameters(unitType, RN_arr) # p0=f(RN), tau1~d0=f(p0)
#                
#                tau1_arr,tau2_arr,KSE_arr,AM_arr,LM_arr,cv_arr,phi1_arr,phi3_arr,\
#                C1i_arr,C1n1_arr,C1n4_arr,C2i_arr,C2n1_arr,C2n4_arr,C3_arr,C4_arr,C5_arr,\
#                alpha_i_arr,beta_arr,gamma_arr,g1_arr,g2_arr,a0_arr,b0_arr,c0_arr,d0_arr = PG.interpParameters(unitType, p0_arr) # p0=f(RN), tau1~d0=f(p0)

        # 241211
        else:
            if(unitType=='musclefibers'):
                MF_arr = np.arange(1, unitCount+1)
                if(p0_file != None): # p0 from file
                    p0_arr = PG.importParameters('p0', p0_file)
                else: # p0 from function
                    p0_arr = PG.genParameters('p0', p0_func, MF_arr, p0_params, unitCount, p0_range)    # f(MF#)
            
            elif(unitType=='motorunit'):
                if(p0_file != None): # p0 from file
                    p0_arr = PG.importParameters('p0', p0_file)
                else: # p0 from function
                    p0_arr = PG.genParameters('p0', p0_func, RN_arr, p0_params, unitCount, p0_range)    # f(RN)

            # default valeus of range parameters
            tau1_arr,tau2_arr,KSE_arr,AM_arr,LM_arr,cv_arr,phi1_arr,phi3_arr,\
            C1i_arr,C1n1_arr,C1n4_arr,C2i_arr,C2n1_arr,C2n4_arr,C3_arr,C4_arr,C5_arr,\
            alpha_i_arr,beta_arr,gamma_arr,g1_arr,g2_arr,a0_arr,b0_arr,c0_arr,d0_arr = PG.interpParameters(unitType, p0_arr) # p0=f(RN), tau1~d0=f(p0)

            # user-defined values of range parameters
            if (tau1_file != None): # tau1 from file
                tau1_arr = PG.importParameters('tau1', tau1_file)
            elif (tau1_func != None): # tau1 from function
                tau1_arr = PG.genParameters('tau1', tau1_func, p0_arr, tau1_params, unitCount, tau1_range) # f(p0)
            
            if (tau2_file != None): # tau2 from file
                tau2_arr = PG.importParameters('tau2', tau2_file)
            elif (tau2_func != None): # tau2 from function
                tau2_arr = PG.genParameters('tau2', tau2_func, p0_arr, tau2_params, unitCount, tau2_range) # f(p0)
            
            if (KSE_file != None): # KSE from file
                KSE_arr = PG.importParameters('KSE', KSE_file)
            elif (KSE_func != None): # KSE from function
                KSE_arr = PG.genParameters('KSE', KSE_func, p0_arr, KSE_params, unitCount, KSE_range) # f(p0)
            
            if (AM_file != None): # AM from file
                AM_arr = PG.importParameters('AM', AM_file)
            elif (AM_func != None): # AM from function
                AM_arr = PG.genParameters('AM', AM_func, p0_arr, AM_params, unitCount, AM_range) # f(p0)
            
            if (LM_file != None): # LM from file
                LM_arr = PG.importParameters('LM', LM_file)
            elif (LM_func != None): # LM from function
                LM_arr = PG.genParameters('LM', LM_func, p0_arr, LM_params, unitCount, LM_range) # f(p0)
            
            if (cv_file != None): # cv from file
                cv_arr = PG.importParameters('cv', cv_file)
            elif (cv_func != None): # cv from function
                cv_arr = PG.genParameters('cv', cv_func, p0_arr, cv_params, unitCount, cv_range) # f(p0)

            if (phi1_file != None): # phi1 from file
                phi1_arr = PG.importParameters('phi1', phi1_file)
            elif (phi1_func != None): # phi1 from function
                phi1_arr = PG.genParameters('phi1', phi1_func, p0_arr, phi1_params, unitCount, phi1_range)  # f(p0)
            
            if (phi3_file != None): # phi3 from file
                phi3_arr = PG.importParameters('phi3', phi3_file)
            elif (phi3_func != None): # phi3 from function
                phi3_arr = PG.genParameters('phi3', phi3_func, p0_arr, phi3_params, unitCount, phi3_range)  # f(p0)

            if (C1i_file != None): # C1i from file
                C1i_arr = PG.importParameters('C1i', C1i_file)
            elif (C1i_func != None): # C1i from function
                C1i_arr = PG.genParameters('C1i', C1i_func, p0_arr, C1i_params, unitCount, C1i_range)  # f(p0)
            
            if (C1n1_file != None): # C1n1 from file
                C1n1_arr = PG.importParameters('C1n1', C1n1_file)
            elif (C1n1_func != None): # C1n1 from function
                C1n1_arr = PG.genParameters('C1n1', C1n1_func, p0_arr, C1n1_params, unitCount, C1n1_range)  # f(p0)
            
            if (C1n4_file != None): # C1n4 from file
                C1n4_arr = PG.importParameters('C1n4', C1n4_file)
            elif (C1n4_func != None): # C1n4 from function
                C1n4_arr = PG.genParameters('C1n4', C1n4_func, p0_arr, C1n4_params, unitCount, C1n4_range)  # f(p0)

            if (C2i_file != None): # C2i from file
                C2i_arr = PG.importParameters('C2i', C2i_file)
            elif (C2i_func != None): # C2i from function
                C2i_arr = PG.genParameters('C2i', C2i_func, p0_arr, C2i_params, unitCount, C2i_range)  # f(p0)
            
            if (C2n1_file != None): # C2n1 from file
                C2n1_arr = PG.importParameters('C2n1', C2n1_file)
            elif (C2n1_func != None): # C2n1 from function
                C2n1_arr = PG.genParameters('C2n1', C2n1_func, p0_arr, C2n1_params, unitCount, C2n1_range)  # f(p0)
            
            if (C2n4_file != None): # C2n4 from file
                C2n4_arr = PG.importParameters('C2n4', C1n4_file)
            elif (C2n4_func != None): # C2n4 from function
                C2n4_arr = PG.genParameters('C2n4', C2n4_func, p0_arr, C2n4_params, unitCount, C2n4_range)  # f(p0)

            if (C3_file != None): # C3 from file
                C3_arr = PG.importParameters('C3', C3_file)
            elif (C3_func != None): # C3 from function
                C3_arr = PG.genParameters('C3', C3_func, p0_arr, C3_params, unitCount, C3_range)  # f(p0)
            
            if (C4_file != None): # C4 from file
                C4_arr = PG.importParameters('C4', C4_file)
            elif (C4_func != None): # C4 from function
                C4_arr = PG.genParameters('C4', C4_func, p0_arr, C4_params, unitCount, C4_range)  # f(p0)
            
            if (C5_file != None): # C5 from file
                C5_arr = PG.importParameters('C5', C5_file)
            elif (C5_func != None): # C5 from function
                C5_arr = PG.genParameters('C5', C5_func, p0_arr, C5_params, unitCount, C5_range)  # f(p0)

            if (alpha_i_file != None): # alpha_i from file
                alpha_i_arr = PG.importParameters('alpha_i', alpha_i_file)
            elif (alpha_i_func != None): # alpha_i from function
                alpha_i_arr = PG.genParameters('alpha_i', alpha_i_func, p0_arr, alpha_i_params, unitCount, alpha_i_range)  # f(p0)
            
            if (beta_file != None): # beta from file
                beta_arr = PG.importParameters('beta', beta_file)
            elif (beta_func != None): # beta from function
                beta_arr = PG.genParameters('beta', beta_func, p0_arr, beta_params, unitCount, beta_range)  # f(p0)
            
            if (gamma_file != None): # gamma from file
                gamma_arr = PG.importParameters('gamma', gamma_file)
            elif (gamma_func != None): # gamma from function
                gamma_arr = PG.genParameters('gamma', gamma_func, p0_arr, gamma_params, unitCount, gamma_range)  # f(p0)

            if (g1_file != None): # g1 from file
                g1_arr = PG.importParameters('g1', beta_file)
            elif (g1_func != None): # g1 from function
                g1_arr = PG.genParameters('g1', g1_func, p0_arr, g1_params, unitCount, g1_range)  # f(p0)
            
            if (g2_file != None): # g2 from file
                g2_arr = PG.importParameters('g2', g2_file)
            elif (g2_func != None): # g2 from function
                g2_arr = PG.genParameters('g2', g2_func, p0_arr, g2_params, unitCount, g2_range)  # f(p0)
            
            if (a0_file != None): # a0 from file
                a0_arr = PG.importParameters('a0', a0_file)
            elif (a0_func != None): # a0 from function
                a0_arr = PG.genParameters('a0', a0_func, p0_arr, a0_params, unitCount, a0_range)  # f(p0)

            if (b0_file != None): # b0 from file
                b0_arr = PG.importParameters('b0', b0_file)
            elif (b0_func != None): # b0 from function
                b0_arr = PG.genParameters('b0', b0_func, p0_arr, b0_params, unitCount, b0_range)  # f(p0)

            if (c0_file != None): # c0 from file
                c0_arr = PG.importParameters('c0', c0_file)
            elif (c0_func != None): # c0 from function
                c0_arr = PG.genParameters('c0', c0_func, p0_arr, c0_params, unitCount, c0_range)  # f(p0)

            if (d0_file != None): # d0 from file
                d0_arr = PG.importParameters('d0', d0_file)
            elif (d0_func != None): # d0 from function
                d0_arr = PG.genParameters('d0', d0_func, p0_arr, d0_params, unitCount, d0_range)  # f(p0)
            # 241211

            # calculartion of phi2, phi4, and g3
            phi2_arr = 1+8*phi1_arr
            phi4_arr = 1+8*phi3_arr
            g3_arr = 1-np.exp(-((-8-g1_arr)/g2_arr)**2)

            num = 1
            for p0,a0,c0,tau1,tau2,KSE,AM,LM,cv,phi1,phi2,phi3,phi4,C1i,C1n1,C1n4,C2i,C2n1,C2n4,C3,C4,C5,alpha_i,beta,gamma,g1,g2,g3,b0,d0 \
                    in zip(p0_arr,a0_arr,c0_arr,tau1_arr,tau2_arr,KSE_arr,AM_arr,LM_arr,cv_arr,phi1_arr,phi2_arr,phi3_arr,phi4_arr,
                           C1i_arr,C1n1_arr,C1n4_arr,C2i_arr,C2n1_arr,C2n4_arr,C3_arr,C4_arr,C5_arr,alpha_i_arr,beta_arr,gamma_arr,
                           g1_arr,g2_arr,g3_arr,b0_arr,d0_arr):
                # params_table.set_value('p0', col, p0) # CHJ
                params_table.at['p0', col] = p0
                # params_table.set_value('a0', col, a0) # CHJ
                params_table.at['a0', col] = a0
                # params_table.set_value('c0', col, c0) # CHJ
                params_table.at['c0', col] = c0
                # params_table.set_value('tau1', col, tau1) # CHJ
                params_table.at['tau1', col] = tau1
                # params_table.set_value('tau2', col, tau2) # CHJ
                params_table.at['tau2', col] = tau2
                # params_table.set_value('Kse', col, KSE) # CHJ
                params_table.at['Kse', col] = KSE
                # params_table.set_value('AM', col, AM) # CHJ
                params_table.at['AM', col] = AM
                # params_table.set_value('LM', col, LM) # CHJ
                params_table.at['LM', col] = LM
                # params_table.set_value('cv', col, cv) # CHJ
                params_table.at['cv', col] = cv
                # params_table.set_value('phi1', col, phi1) # CHJ
                params_table.at['phi1', col] = phi1
                # params_table.set_value('phi2', col, phi2) # CHJ
                params_table.at['phi2', col] = phi2
                # params_table.set_value('phi3', col, phi3) # CHJ
                params_table.at['phi3', col] = phi3
                # params_table.set_value('phi4', col, phi4) # CHJ
                params_table.at['phi4', col] = phi4
                # params_table.set_value('C1i', col, C1i) # CHJ
                params_table.at['C1i', col] = C1i
                # params_table.set_value('C1n1', col, C1n1) # CHJ
                params_table.at['C1n1', col] = C1n1
                # params_table.set_value('C1n4', col, C1n4) # CHJ
                params_table.at['C1n4', col] = C1n4
                # params_table.set_value('C2i', col, C2i) # CHJ
                params_table.at['C2i', col] = C2i
                # params_table.set_value('C2n1', col, C2n1) # CHJ
                params_table.at['C2n1', col] = C2n1
                # params_table.set_value('C2n4', col, C2n4) # CHJ
                params_table.at['C2n4', col] = C2n4
                # params_table.set_value('C3', col, C3) # CHJ
                params_table.at['C3', col] = C3
                # params_table.set_value('C4', col, C4) # CHJ
                params_table.at['C4', col] = C4
                # params_table.set_value('C5', col, C5) # CHJ
                params_table.at['C5', col] = C5
                # params_table.set_value('alpha_i', col, alpha_i) # CHJ
                params_table.at['alpha_i', col] = alpha_i
                # params_table.set_value('beta', col, beta) # CHJ
                params_table.at['beta', col] = beta
                # params_table.set_value('gamma', col, gamma) # CHJ
                params_table.at['gamma', col] = gamma
                # params_table.set_value('g1', col, g1) # CHJ
                params_table.at['g1', col] = g1
                # params_table.set_value('g2', col, g2) # CHJ
                params_table.at['g2', col] = g2
                # params_table.set_value('g3', col, g3) # CHJ
                params_table.at['g3', col] = g3
                # params_table.set_value('b0', col, b0) # CHJ
                params_table.at['b0', col] = b0
                # params_table.set_value('d0', col, d0) # CHJ
                params_table.at['d0', col] = d0
                params = params_table[:][col].values
                mf_params.append(params.copy())  # ndarray
                
                # display of range parameter values
                # print 'MTU %d : p0 = %.3f, tau1 = %.3f, tau2 = %.3f, KSE = %.3f, AM = %.3f, LM = %.3f, cv = %.3f, ' \
                #       'phi1 = %.4f, phi2 = %.3f, phi3 = %.4f, phi4 = %.3f, ' \
                #       'C1i = %.3f, C1n1 = %.3f, C1n4 = %.4f, C2i = %.3f, C2n1 = %.3f, C2n4 = %.4f, ' \
                #       'C3 = %.3f, C4 = %.3f, C5 = %.3f, ' \
                #       'alpha_i = %.3f, beta = %.3f, gamma = %.4f, ' \
                #       'g1 = %.3f, g2 = %.3f, g3 = %.4f, a0 = %.3f, b0 = %.3f, c0 = %.3f, d0 = %.3f' \
                #       % (num, p0, tau1, tau2, KSE, AM, LM, cv, phi1, phi2, phi3, phi4, C1i, C1n1, C1n4, C2i, C2n1, C2n4, C3, C4, C5, alpha_i, beta, gamma, g1, g2, g3, a0, b0, c0, d0) # CHJ
                print('MTU %d : p0 = %.3f, tau1 = %.3f, tau2 = %.3f, KSE = %.3f, AM = %.3f, LM = %.3f, cv = %.3f, ' \
                      'phi1 = %.4f, phi2 = %.3f, phi3 = %.4f, phi4 = %.3f, ' \
                      'C1i = %.3f, C1n1 = %.3f, C1n4 = %.4f, C2i = %.3f, C2n1 = %.3f, C2n4 = %.4f, ' \
                      'C3 = %.3f, C4 = %.3f, C5 = %.3f, ' \
                      'alpha_i = %.3f, beta = %.3f, gamma = %.4f, ' \
                      'g1 = %.3f, g2 = %.3f, g3 = %.4f, a0 = %.3f, b0 = %.3f, c0 = %.3f, d0 = %.3f' \
                      % (num, p0, tau1, tau2, KSE, AM, LM, cv, phi1, phi2, phi3, phi4, C1i, C1n1, C1n4, C2i, C2n1, C2n4, C3, C4, C5, alpha_i, beta, gamma, g1, g2, g3, a0, b0, c0, d0))
                num += 1

        if(unitType=='musclefibers'):
            MP.setParameters(mf_params)
        elif(unitType=='motorunit'):
            MP.setParameters(mn_params, mf_params)


## Range parameter display
def plot_params(param, unitType=None):

    global MP, PG
    unitCount=MP.getUnitCount()
    MP_arr = np.arange(1, unitCount+1)
    PG.plotParameters(MP_arr, param, unitType)


## Simulation time setting
def setSimulTimes(t_start=0, t_stop=10000, t_dt=0.1):
    
    global MP
    MP.setSimulTimes(t_start, t_stop, t_dt)  


## State variable initializaiton
def setInitialValues(ivalue_1, ivalue_2=None):
       
    global MP
    unitType=MP.getUnitType()

    if(unitType=='motorunit'):
        MP.setInitialValues(ivalue_1, ivalue_2)
        
    else:
        MP.setInitialValues(ivalue_1)
 

## Intracellular current from built-in function for motoneuron  
def genNeuronInputSignals(signalType, period=5000., amplitude=20., offset=0., ivalue=0., heav_param=Is_heav_param):
    """
        Is = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
           + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
           + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
           + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
           + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
    """
    
    global MP, NSG
    t_start, t_stop, t_dt =MP.getSimulTimes()
    NSG.genSignal(signalType, period, amplitude, offset, ivalue, t_stop, t_dt, heav_param)
    

## Intracellular current from data file for motoneuron 
def importNeuronInputSignals(filePath):
    global MP, NSG
    NSG.importSignalFile(filePath)
    

## Intracellular current assignment for motoneuron population
def setNeuronInputSignals(MN_range_1, MN_range_2):
    global MP, NSG
    NSG.setSignal(MP, MN_range_1, MN_range_2)


## Synaptic input from built-in function for motoneuron 
def genSynConSignals(department, synType, signalType='Ramp', heav_param=Syn_heav_param, iValue=0., pValue=0.1, period=10000, tau=0., std_max=0., noise=False):
    """
        SynCon = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
           + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
           + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
           + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
           + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
    """
    global MP, s_SCSG, d_SCSG
    t_start, t_stop, t_dt =MP.getSimulTimes()
        
    if(department == 'Soma'):
        s_SCSG.genSignal(synType, signalType, heav_param, iValue, pValue, period, tau, std_max, noise, t_stop, t_dt)

    elif(department == 'Dendrite'):
        d_SCSG.genSignal(synType, signalType, heav_param, iValue, pValue, period, tau, std_max, noise, t_stop, t_dt)
    

## Synaptic input from data file for motoneuron   
def importSynConSignals(department, synType, filePath):
    global MP, s_SCSG, d_SCSG

    if(department == 'Soma'):
        s_SCSG.importSignalFile(synType, filePath)

    if(department == 'Dendrite'):
        d_SCSG.importSignalFile(synType, filePath)


## Synaptic input assignment for motoneuron population
def setSynConSignals(department, MN_range_1, MN_range_2):
    global MP, s_SCSG, d_SCSG
    
    if(department == 'Soma'):
        s_SCSG.setSignal(MP, MN_range_1, MN_range_2, department)

    elif(department == 'Dendrite'):
        d_SCSG.setSignal(MP, MN_range_1, MN_range_2, department)


## Impulse current input from built-in function for muscle-tendon unit
def genSpikeSignals(signalType, t_start=100., t_stop=1400., freq=100., scale=0.):
    global MP, SSG
    t_1, t_2, t_dt =MP.getSimulTimes()
    SSG.genSignal(signalType, t_start, t_stop, freq, scale, t_dt)
    

## Impulse current input from data file for muscle-tendon unit        
def importSpikeSignals(filePath, t_dt=0.1):
    global MP, SSG
    SSG.importSignalFile(filePath, t_dt)
        

## Impulse current input assignment for muscle-tendon unit population   
def setSpikeSignals(MF_range_1, MF_range_2):
    global MP, SSG
    SSG.setSignal(MP, MF_range_1, MF_range_2)


## Length variation from built-in function for muscle-tendon unit 
def genMuscleLengthSignals(signalType, itime=-8., ftime=0., ivalue=0., fvalue=0.):
    global MP, MSG
    t_start, t_stop, t_dt = MP.getSimulTimes()
    MSG.genSignal(signalType, itime, ftime, ivalue, fvalue, t_stop, t_dt)
    

## Length variation from data file for muscle-tendon unit
def importMuscleLengthSignals(filePath):
    global MP, MSG
    MSG.importSignalFile(filePath)
    

## Length variation assignment for muscle-tendon unit population    
def setMuscleLengthSignals(MF_range_1, MF_range_2):
    global MP, MSG
    MSG.setSignal(MP, MF_range_1, MF_range_2)


## Intracellular current preview for motoneuron    
def plotNeuronInputSignal():
    global NSG
    NSG.plotSignal()
    

## Synaptic input preview for motoneuron    
def plotSynConSignal(department):
    global s_SCSG, d_SCSG
    if(department == 'Soma'):
        s_SCSG.plotSignal(department)
    elif(department == 'Dendrite'):
        d_SCSG.plotSignal(department)


## Current impulse preview for muscle-tendon unit
def plotSpikeSignal():
    global SSG
    SSG.plotSignal()


## Length preview for muscle-tendon unit
def plotMuscleLengthSignal():
    global MSG
    MSG.plotSignal()


## Computational node selection
def setComputeNode(node_list):
    global PM
    PM=sim.ParallelManager()
    PM.selectComputeNode(node_list) 


## Parallel simulation control    
def runSimulation(num, node_list):
    global MP, PM
    PM.selectJobServer(num, node_list)
    PM.runSimulation(MP)


## Simulation result saving
def saveSimulationResults(savePath, fileName):
    global MP, PM
    PM.saveResultData(MP, savePath, fileName)


## Online simulation result display
def plotSimulResult(scope_list):
    global PM, MP
    PM.plotResultData(MP, scope_list)

    
## Offline simulation result display    
def plotImportData(dirPath, fileName, display, scope_list):
    global MP
    unitCount=MP.getUnitCount()
    unitType=MP.getUnitType()
    
    if(display == 'Individual'):
        nrows=unitCount
        ncols=1
        gs0 = gridspec.GridSpec(nrows, ncols)
    
        for scope in scope_list: # scope_list = vs, vd
    
            for i in range(1, unitCount+1):
                filePath = dirPath + fileName + str('_%03d' %i) + '.csv'
                import_data = pd.read_csv(filePath) # reading from DB. dataFrame

                # 240906
                #x_data=import_data['Time']
                x_data = import_data['Time']/1000
                y_data=import_data[scope]
                ax = plt.subplot(gs0[i-1])
    
                if(scope=='Firing_rate'):
                    x_data=x_data[import_data['Firing_rate'] > 0]
                    y_data=import_data['Firing_rate'][import_data['Firing_rate'] > 0]
                    plt.plot(x_data, y_data, '.')
                else:
                    plt.plot(x_data, y_data)

                ax.set_xlim(0, np.asarray(import_data['Time'])[-1]/1000)
                ax.set_autoscaley_on(True)
                ax.autoscale_view()

                if (unitType == 'motoneuron'):
                    plt.ylabel('MN' + str(i))
                elif (unitType == 'musclefibers'):
                    plt.ylabel('MTU' + str(i))
                elif (unitType == 'motorunit'):
                    plt.ylabel('MU' + str(i))

                plt.grid('on')
    
            plt.xlabel('Time (sec)')
            figure = plt.gcf()
            # figure.canvas.set_window_title(scope) # CHJ
            figure.canvas.manager.set_window_title(scope)
            # plt.show() # CHJ
            plt.show(block=True)
        
    elif(display == 'Combined'):
        for scope in scope_list: # scope_list = vs, vd
            ax = plt.subplot(111)
            ax.set_autoscaley_on(True)
            ax.autoscale_view()
            plt.grid('on')
            
            for i in range(1, unitCount+1):
                filePath = dirPath + fileName + str('_%03d' %i) + '.csv'
                import_data = pd.read_csv(filePath) # reading from DB. dataFrame

                x_data = import_data['Time']/1000
                y_data=import_data[scope]
                
                if(scope=='Firing_rate'):
                    x_data=x_data[import_data['Firing_rate'] > 0]
                    y_data=import_data['Firing_rate'][import_data['Firing_rate'] > 0]

                    if (unitType == 'motoneuron'):
                        plt.plot(x_data, y_data, '.', label='MN' + str(i))
                    elif (unitType == 'musclefibers'):
                        plt.plot(x_data, y_data, '.', label='MTU' + str(i))
                    elif (unitType == 'motorunit'):
                        plt.plot(x_data, y_data, '.', label='MU' + str(i))
                    
                else:
                    if (unitType == 'motoneuron'):
                        plt.plot(x_data, y_data, label='MN' + str(i))
                    elif (unitType == 'musclefibers'):
                        plt.plot(x_data, y_data, label='MTU' + str(i))
                    elif (unitType == 'motorunit'):
                        plt.plot(x_data, y_data, label='MU' + str(i))
                    
            ax.set_xlim(0, np.asarray(import_data['Time'])[-1]/1000)
            plt.xlabel('Time (sec)')
            ax.legend(loc='best')
            figure = plt.gcf()
            # figure.canvas.set_window_title(scope) # CHJ
            figure.canvas.manager.set_window_title(scope)
            # plt.show() # CHJ
            plt.show(block=True)
            
    elif(display == 'Sum'):
        for scope in scope_list: # scope_list = vs, vd
            ax = plt.subplot(111)
            ax.set_autoscaley_on(True)
            ax.autoscale_view()
            plt.grid('on')
            y_sum = np.zeros(1)
            
            for i in range(1, unitCount+1):
                filePath = dirPath + fileName + str('_%03d' %i) + '.csv'
                import_data = pd.read_csv(filePath) # reading from DB. dataFrame
        
                x_data = import_data['Time']/1000
                y_data = import_data[scope]
                y_sum = y_sum + y_data
                
            if(scope=='Firing_rate'):
                x_data = x_data[y_sum > 0]
                y_sum = y_sum[y_sum > 0]
                plt.plot(x_data, y_sum, '.')
            else:
                plt.plot(x_data, y_sum)

            ax.set_xlim(0, np.asarray(import_data['Time'])[-1]/1000)
            plt.xlabel('Time (sec)')
            figure = plt.gcf()
            # figure.canvas.set_window_title(scope+' Sum up') # CHJ
            figure.canvas.manager.set_window_title(scope+' Sum up')
            
            # plt.show() # CHJ
            plt.show(block=True)