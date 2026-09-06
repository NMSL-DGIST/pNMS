# -*- coding: utf-8 -*-
# pMUPS (ver3.0)

## Backend for HPC
import matplotlib
# matplotlib.use('qt4agg', warn=False) # CHJ
matplotlib.use('Qt5Agg')
## Backend for HPC

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import matplotlib as mpl # CHJ
from cycler import cycler # CHJ
mpl.rcParams['axes.prop_cycle'] = cycler('color', ['b', 'g', 'r', 'c', 'm', 'y', 'k']) # CHJ

from math import log, exp, pi, asin, sin, tanh, cosh
from scipy import integrate, signal
from scipy.interpolate import splev, splrep, interp2d
import pandas as pd
import platform
import pp


## Parent class for motoneuron and muscle-tenton unit
class Cell(object):

    def __init__(self, *args, **kwargs):
        
        self.setID(kwargs['uniqueNumber'])
        # print "Cell instance constructed. ID: "+str(self.uniqueNumber) # CHJ
        print("Cell instance constructed. ID: "+str(self.uniqueNumber))
        
    def __del__(self):
        # print "Cell instance destructed. ID: "+str(self.uniqueNumber) # CHJ
        print("Cell instance destructed. ID: "+str(self.uniqueNumber))

    def setCellType(self, cellType):
        self.cellType=cellType

    def getCellType(self):
        return self.cellType

    def setID(self, num):
        self.uniqueNumber=num

    def getID(self):
        return self.uniqueNumber

    def getHostname(self):
        self.hostname = platform.node()
        return self.hostname

    def setSimulTime(self, t_start, t_stop, t_dt):
        self.t_start=t_start
        self.t_stop=t_stop
        self.t_dt=t_dt
        
        # print "Simulation time changed. ID: "+str(self.getID()) # CHJ
        print("Simulation time changed. ID: "+str(self.getID()))
        # print "start time: "+str(self.t_start)+" (msec)" # CHJ
        print("start time: "+str(self.t_start)+" (msec)")
        # print "stop time: "+str(self.t_stop)+" (msec)" # CHJ
        print("stop time: "+str(self.t_stop)+" (msec)")
        # print "time interval: "+str(self.t_dt)+" (msec) # CHJ
        print("time interval: "+str(self.t_dt)+" (msec)")
        

    def getSimulTime(self):
        return self.t_stop, self.t_dt

    def setParameter(self, parameter):
        self.parameter = parameter

    def getParameter(self):
        return self.parameter

    def setInitialValue(self, ivalue):
        self.ivalue = ivalue
        # print str(self.getCellType())+" initial values set." # CHJ
        print(str(self.getCellType())+" initial values set.")

    def getInitialValue(self):
        return self.ivalue
        

## Child class of Cell class for motoneurons
class Motoneuron(Cell, object):
    def __init__(self, *args, **kwargs):
        super(Motoneuron, self).__init__(*args, **kwargs)
        self.setCellType('motoneuron')
        
        # print "Motoneuron instance constructed. ID: "+str(self.uniqueNumber) # CHJ
        print("Motoneuron instance constructed. ID: "+str(self.uniqueNumber))

    def __del__(self):
        # print "Motoneuron instance destructed. ID: "+str(self.uniqueNumber) # CHJ
        print("Motoneuron instance destructed. ID: "+str(self.uniqueNumber))

    def heav(self, x):
        return np.where(x < 0, 0, 1)

    def model(self, t, y):
        # state variables
        vs, sca, snam, snah, snapm, skdr, scam, scah, shm, vd,  dca, dcal, dnam, dnah, dnapm, dkdr, dkcam, dcam, dcah, dhm = y
        
        # model parameters
        gms,gmd,gc,cmd,cms, \
        svl,dvl, \
        sf,skca,salpha,sCAo,sEca, \
        sgna,svna,sanamc,sanamv,sanama,sanamb,sbnamc,sbnamv,sbnama,sbnamb,snahth,snahslp,snahv,snaha,snahb,snahc, \
        sgnap,svna,sanapmc,sanapmv,sanapma,sanapmb,sbnapmc,sbnapmv,sbnapma,sbnapmb, \
        sgkdr,svk,skdrth,skdrslp,skdrv,skdra,skdrb,skdrc, \
        sgkca,svk,skd, \
        sgca,scamth,scamslp,scamtau,scahth,scahslp,scahtau, \
        sgh,svh,shth,shslp,shtau, \
        svesyn, \
        svisyn, \
        df,dkca,dalpha,dCAo,dEca, \
        dgcal,dcalth,dcalslp,dcaltau, \
        dgna,dvna,danamc,danamv,danama,danamb,dbnamc,dbnamv,dbnama,dbnamb,dnahth,dnahslp,dnahv,dnaha,dnahb,dnahc, \
        dgnap,dvna,danapmc,danapmv,danapma,danapmb,dbnapmc,dbnapmv,dbnapma,dbnapmb, \
        dgkdr,dvk,dkdrth,dkdrslp,dkdrv,dkdra,dkdrb,dkdrc, \
        dgkca,dvk,dkcamth,dkcamslp,dkcamtau, \
        dgca,dcamth,dcamslp,dcamtau,dcahth,dcahslp,dcahtau, \
        dgh,dvh,dhth,dhslp,dhtau, \
        dvesyn, \
        dvisyn=self.parameter
        
        # constants
        R=8.31441
        Temp=309.15
        Zca=2
        Fe=96485.309
        const = 1000*R*Temp/Zca/Fe
        p = 0.492
            
        # intracellular current stimulation
        if(self.is_sigType=='Step'):
            i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.is_heav_param
            I_s = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                         + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                         + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                         + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                         + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
        elif(self.is_sigType=='Ramp'):
            I_s=self.is_amp-((self.is_amp-self.is_iv)/self.is_period)*abs(t-self.is_period)+self.is_offset
        elif(self.is_sigType=='Import'):
            I_s=np.interp(t, self.is_time,  self.Is, 0, 0)
        elif(self.is_sigType=='sine'):
            I_s=self.is_amp*sin(2*pi*t/self.is_period)+self.is_offset
        elif(self.is_sigType=='square'):
            x_t=t%self.is_period
            if(x_t<(self.is_period/2)):
                I_s=1*self.is_amp+self.is_offset
            elif(x_t>=self.is_period/2):
                I_s=0+self.is_offset
        I_s = I_s/3.1576 # unit conversion from nA to mA/cm^2

        # synaptic conductance variation
        sgesyn = np.interp(t, self.sesyn_t, self.Se_syncon, 0, 0)
        sgisyn = np.interp(t, self.sisyn_t, self.Si_syncon, 0, 0)
        dgesyn = np.interp(t, self.desyn_t, self.De_syncon, 0, 0)
        dgisyn = np.interp(t, self.disyn_t, self.Di_syncon, 0, 0)
                        
        # active currents at the somatic compartment
        # intracellular calcium dynamics
        if(sca < 1.e-100):
            sca = 1.e-100
        # Ca2+ reversal potential
        if(sEca != 0.): # constant
            svca = sEca
        else: # dynamic
            svca = const*log(sCAo/sca)-70
        
        # N-Type Ca2+ current
        sica=sgca*scam**2*scah*(vs-svca)     
        # fast Na+ current
        sina=sgna*snam**3*snah*(vs-svna)
        # persistent Na+ current
        sinap=sgnap*snapm**3*(vs-svna)
        # HCN current
        sih=sgh*shm*(vs-svh)
        # synaptic current
        siesyn=sgesyn*(vs-svesyn)
        siisyn=sgisyn*(vs-svisyn)
        sisyn=siesyn+siisyn
        # total inward current 
        ins=-sina-sica-sinap-sih-sisyn

        # delayed rectifier K+ current
        sikdr=sgkdr*skdr**4*(vs-svk)
        # Ca2+ dependent K+ current
        sikca=sgkca*(sca/(sca+skd))*(vs-svk)
        # total outward current
        outs=-sikdr-sikca
        
        # active currents at the dendritic compartment
        # intracellular calcium dynamics
        if(dca < 1.e-100):
            dca = 1.e-100        
        # Ca2+ reversal potential
        if(dEca != 0.): # constant
            dvca = dEca
        else: # dynamic
            dvca=const*log(dCAo/dca)-70

        # L-Type Ca2+ current
        dical=dgcal*dcal*(vd-dvca)
        # N-Type Ca2+ current
        dica=dgca*dcam**2*dcah*(vd-dvca)
        # fast Na+ current
        dina=dgna*dnam**3*dnah*(vd-dvna)
        # persistent Na+ current
        dinap=dgnap*dnapm**3*(vd-dvna)
        # HCN current
        dih=dgh*dhm*(vd-dvh)
        # synaptic current 
        diesyn=dgesyn*(vd-dvesyn)
        diisyn=dgisyn*(vd-dvisyn)
        disyn=diesyn+diisyn
        # total inward current
        ind=-dical-dina-dica-dinap-dih-disyn

        # delayed rectifier K+ current
        dikdr=dgkdr*dkdr**4*(vd-dvk)
        # Ca2+ dependent K+ current
        dikca=dgkca*dkcam*(vd-dvk)
        # total outward current
        outd=-dikdr-dikca

        # model equations
        n = len(y)      
        # dydt=range(n) # CHJ
        dydt = list(range(n))

        # d(vs)/dt
        dydt[0]=(I_s+ins+outs-gms*(vs-svl)+gc*(vd-vs)/p)/cms
        
        # d(sca)/dt
        dydt[1]=-sf*salpha*sica-sf*skca*sca
        
        # d(snam)/dt
        alpha_snafm=sanamc*(vs-sanamv)/(exp(-(vs-sanamv)/sanama)+sanamb)
        beta_snafm=sbnamc*(vs-sbnamv)/(exp((vs-sbnamv)/sbnama)+sbnamb)
        dydt[2]=alpha_snafm*(1-snam)-beta_snafm*snam

        # d(snah)/dt
        snahinf=1.0/(1.0+exp((vs-snahth)/snahslp))
        snahtau=snahc/(exp((vs-snahv)/snaha)+exp(-(vs-snahv)/snahb))
        dydt[3]=(snahinf-snah)/snahtau

        # d(snapm))/dt
        alpha_snapm=sanapmc*(vs-sanapmv)/(exp(-(vs-sanapmv)/sanapma)+sanapmb)
        beta_snapm=sbnapmc*(vs-sbnapmv)/(exp((vs-sbnapmv)/sbnapma)+sbnapmb)
        dydt[4]=alpha_snapm*(1-snapm)-beta_snapm*snapm
    
        # d(skdr)/dt
        skdrinf=1.0/(1.0+exp(-(vs-skdrth)/skdrslp))
        skdrtau=skdrc/(exp((vs-skdrv)/skdra)+exp(-(vs-skdrv)/skdrb))
        dydt[5]=(skdrinf-skdr)/skdrtau

        # d(scam)/dt
        scaminf=1.0/(1.0+exp(-(vs-scamth)/scamslp))
        dydt[6]=(scaminf-scam)/scamtau

        # d(scah)/dt
        scahinf=1.0/(1.0+exp((vs-scahth)/scahslp))
        dydt[7]=(scahinf-scah)/scahtau
        
        # d(shm)/dt 
        shminf=1/(1+exp((vs+shth)/shslp))
        dydt[8]=(shminf-shm)/shtau


        # d(vd)/dt 
        dydt[9]=(ind+outd-gmd*(vd-dvl)+gc*(vs-vd)/(1.0-p))/cmd
        
        # total calcium current
        sdica = dical+dica
        # d(dca)/dt 
        dydt[10]=-df*dalpha*sdica-df*dkca*dca

        # d(dcal)/dt
        dcalinf=1.0/(1.0+exp(-(vd-dcalth)/dcalslp))
        dydt[11]=(dcalinf-dcal)/dcaltau
        
        # d(dnam)/dt 
        alpha_dnafm=danamc*(vd-danamv)/(exp(-(vd-danamv)/danama)+danamb)
        beta_dnafm=dbnamc*(vd-dbnamv)/(exp((vd-dbnamv)/dbnama)+dbnamb)
        dydt[12]=alpha_dnafm*(1-dnam)-beta_dnafm*dnam 
    
        # d(dnah)/dt 
        dnahinf=1.0/(1.0+exp((vd-dnahth)/dnahslp))
        dnahtau=dnahc/(exp((vd-dnahv)/dnaha)+exp(-(vd-dnahv)/dnahb))
        dydt[13]=(dnahinf-dnah)/dnahtau
        
        # d(dnapm)/dt 
        alpha_dnapm=danapmc*(vd-danapmv)/(exp(-(vd-danapmv)/danapma)+danapmb)
        beta_dnapm=dbnapmc*(vd-dbnapmv)/(exp((vd-dbnapmv)/dbnapma)+dbnapmb)
        dydt[14]=alpha_dnapm*(1-dnapm)-beta_dnapm*dnapm
        
        # d(dkdr)/dt 
        dkdrinf=1.0/(1.0+exp(-(vd-dkdrth)/dkdrslp))
        dkdrtau=dkdrc/(exp((vd-dkdrv)/dkdra)+exp(-(vd-dkdrv)/dkdrb))
        dydt[15]=(dkdrinf-dkdr)/dkdrtau
        
        # d(dkcam)/dt
        dkcaminf=1/(1+(dkcamth/dca)**dkcamslp)
        dydt[16]=(dkcaminf-dkcam)/dkcamtau

        # d(dcam)/dt 
        dcaminf=1.0/(1.0+exp(-(vd-dcamth)/dcamslp))
        dydt[17]=(dcaminf-dcam)/dcamtau

        # d(dcah)/dt 
        dcahinf=1.0/(1.0+exp((vd-dcahth)/dcahslp))
        dydt[18]=(dcahinf-dcah)/dcahtau

        # d(dhm)/dt 
        dhminf=1/(1+exp((vd+dhth)/dhslp))
        dydt[19]=(dhminf-dhm)/dhtau
    
        return dydt
            
    def solModel(self):
        # model parameters
        gms,gmd,gc,cmd,cms, \
        svl,dvl, \
        sf,skca,salpha,sCAo,sEca, \
        sgna,svna,sanamc,sanamv,sanama,sanamb,sbnamc,sbnamv,sbnama,sbnamb,snahth,snahslp,snahv,snaha,snahb,snahc, \
        sgnap,svna,sanapmc,sanapmv,sanapma,sanapmb,sbnapmc,sbnapmv,sbnapma,sbnapmb, \
        sgkdr,svk,skdrth,skdrslp,skdrv,skdra,skdrb,skdrc, \
        sgkca,svk,skd, \
        sgca,scamth,scamslp,scamtau,scahth,scahslp,scahtau, \
        sgh,svh,shth,shslp,shtau, \
        svesyn, \
        svisyn, \
        df,dkca,dalpha,dCAo,dEca, \
        dgcal,dcalth,dcalslp,dcaltau, \
        dgna,dvna,danamc,danamv,danama,danamb,dbnamc,dbnamv,dbnama,dbnamb,dnahth,dnahslp,dnahv,dnaha,dnahb,dnahc, \
        dgnap,dvna,danapmc,danapmv,danapma,danapmb,dbnapmc,dbnapmv,dbnapma,dbnapmb, \
        dgkdr,dvk,dkdrth,dkdrslp,dkdrv,dkdra,dkdrb,dkdrc, \
        dgkca,dvk,dkcamth,dkcamslp,dkcamtau, \
        dgca,dcamth,dcamslp,dcamtau,dcahth,dcahslp,dcahtau, \
        dgh,dvh,dhth,dhslp,dhtau, \
        dvesyn, \
        dvisyn=self.parameter

        # constants
        R=8.31441
        Temp=309.15
        Zca=2
        Fe=96485.309
        const = 1000*R*Temp/Zca/Fe
        
        SpikeTimes=[]
        FiringRate=[]
    
        # numerical solver setting
        r1= integrate.ode(self.model).set_integrator('lsoda', max_step=10*self.t_dt)
        r1.set_initial_value(self.ivalue, self.t_start) 
       
        # arrays for integration results
        num_steps = int(np.floor((self.t_stop - self.t_start)/self.t_dt) + 1)
        FR = np.zeros(num_steps)
        T = np.zeros(num_steps)
        IS = np.zeros(num_steps)
        VS = np.zeros(num_steps)
        SCA = np.zeros(num_steps)
        SVCA = np.zeros(num_steps)
        SNAI = np.zeros(num_steps)
        SNAM = np.zeros(num_steps)
        SNAH = np.zeros(num_steps)
        SNAPI = np.zeros(num_steps)
        SNAPM = np.zeros(num_steps) 
        SKDRI = np.zeros(num_steps)
        SKDR = np.zeros(num_steps)
        SKCAI = np.zeros(num_steps)
        SCAI = np.zeros(num_steps)
        SCAM = np.zeros(num_steps) 
        SCAH = np.zeros(num_steps)
        SHI = np.zeros(num_steps)
        SHM = np.zeros(num_steps)
        VD = np.zeros(num_steps) 
        DCA = np.zeros(num_steps)
        DVCA = np.zeros(num_steps)
        DCALI = np.zeros(num_steps)
        DCAL = np.zeros(num_steps) 
        DNAI = np.zeros(num_steps)
        DNAM = np.zeros(num_steps)
        DNAH = np.zeros(num_steps)
        DNAPI = np.zeros(num_steps)
        DNAPM = np.zeros(num_steps)
        DKDRI = np.zeros(num_steps)
        DKDR = np.zeros(num_steps)
        DKCAI = np.zeros(num_steps)
        DKCAM = np.zeros(num_steps)
        DCAI = np.zeros(num_steps)
        DCAM = np.zeros(num_steps)
        DCAH = np.zeros(num_steps)
        DHI = np.zeros(num_steps)
        DHM = np.zeros(num_steps)
        SESYNI = np.zeros(num_steps)
        SGESYN = np.zeros(num_steps)
        SISYNI = np.zeros(num_steps)
        SGISYN = np.zeros(num_steps)
        DESYNI = np.zeros(num_steps)
        DGESYN = np.zeros(num_steps)
        DISYNI = np.zeros(num_steps)
        DGISYN = np.zeros(num_steps)
             
        # initialization of state variables 
        IS[0]= self.Is[0]
        VS[0]=self.ivalue[0]
        SCA[0]=self.ivalue[1]
        if(sEca != 0.):
            SVCA[0] = sEca
        else:
            SVCA[0]=(const*log(sCAo/SCA[0]))-70
        SNAM[0]=self.ivalue[2]
        SNAH[0]=self.ivalue[3]
        SNAI[0] = sgna*SNAM[0]**3*SNAH[0]*(VS[0]-svna)
        SNAPM[0]=self.ivalue[4]
        SNAPI[0] = sgnap*SNAPM[0]**3*(VS[0]-svna)
        SKDR[0]=self.ivalue[5]
        SKDRI[0] = sgkdr*SKDR[0]**4*(VS[0]-svk)
        SKCAI[0] = sgkca*(SCA[0]/(SCA[0]+skd))*(VS[0]-svk)
        SCAM[0]=self.ivalue[6]
        SCAH[0]=self.ivalue[7]
        SCAI[0] = sgca*SCAM[0]**2*SCAH[0]*(VS[0]-SVCA[0])
        SHM[0]=self.ivalue[8]
        SHI[0] = sgh*SHM[0]*(VS[0]-svh)
        VD[0]=self.ivalue[9]
        DCA[0]=self.ivalue[10]
        if(dEca != 0.):
            DVCA[0] = dEca
        else:
            DVCA[0]=(const*log(dCAo/DCA[0]))-70
        DCAL[0]=self.ivalue[11]
        DCALI[0] = dgcal*DCAL[0]*(VD[0]-DVCA[0])
        DNAM[0]=self.ivalue[12]
        DNAH[0]=self.ivalue[13]
        DNAI[0] = dgna*DNAM[0]**3*DNAH[0]*(VD[0]-dvna)
        DNAPM[0]=self.ivalue[14]
        DNAPI[0] = dgnap*DNAPM[0]**3*(VD[0]-dvna)
        DKDR[0]=self.ivalue[15]
        DKDRI[0] = dgkdr*DKDR[0]**4*(VD[0]-dvk)
        DKCAM[0] = self.ivalue[16]
        DKCAI[0] = dgkca*DKCAM[0]*(VD[0]-dvk)
        DCAM[0]=self.ivalue[17]
        DCAH[0]=self.ivalue[18]
        DCAI[0] = dgca*DCAM[0]**2*DCAH[0]*(VD[0]-DVCA[0])
        DHM[0]=self.ivalue[19]
        DHI[0] = dgh*DHM[0]*(VD[0]-dvh)
        SGESYN[0] = self.Se_syncon[0]
        SESYNI[0] = SGESYN[0]*(VS[0]-svesyn)
        SGISYN[0] = self.Si_syncon[0]
        SISYNI[0] = SGISYN[0]*(VS[0]-svisyn)
        DGESYN[0] = self.De_syncon[0]
        DESYNI[0] = DGESYN[0]*(VD[0]-dvesyn)
        DGISYN[0] = self.Di_syncon[0]
        DISYNI[0] = DGISYN[0]*(VD[0]-dvisyn)
        
        # array for variables that are saved in result files
        resultArray={'Time': T,
                      'Is' : IS,
                      'Firing_rate' : FR,
                      'V_soma': VS,
                      '[Ca]_soma': SCA,
                      'E_Ca_soma': SVCA,
                      'I_Naf_soma': SNAI,
                      'm_Naf_soma': SNAM,
                      'h_Naf_soma': SNAH,
                      'I_Nap_soma': SNAPI,
                      'm_Nap_soma': SNAPM,
                      'I_Kdr_soma': SKDRI,
                      'n_Kdr_soma' : SKDR,
                      'I_Kca_soma': SKCAI,
                      'I_Can_soma': SCAI,
                      'm_Can_soma': SCAM,
                      'h_Can_soma': SCAH,
                      'I_H_soma': SHI,
                      'm_H_soma' : SHM,    
                      'V_dend': VD,
                      '[Ca]_dend': DCA,
                      'E_Ca_dend': DVCA,
                      'I_Cal_dend': DCALI,
                      'm_Cal_dend': DCAL,
                      'I_Naf_dend': DNAI,
                      'm_Naf_dend': DNAM,
                      'h_Naf_dend': DNAH,
                      'I_Nap_dend': DNAPI,
                      'm_Nap_dend': DNAPM,
                      'I_Kdr_dend': DKDRI,
                      'n_Kdr_dend' : DKDR,
                      'I_Kca_dend': DKCAI,
                      'm_Kca_dend': DKCAM,
                      'I_Can_dend': DCAI,
                      'm_Can_dend': DCAM,
                      'h_Can_dend': DCAH,
                      'I_H_dend': DHI,
                      'm_H_dend' : DHM,
                      'I_esyn_soma': SESYNI,
                      'G_esyn_soma': SGESYN,
                      'I_isyn_soma': SISYNI,
                      'G_isyn_soma': SGISYN,
                      'I_esyn_dend': DESYNI,
                      'G_esyn_dend': DGESYN, 
                      'I_isyn_dend': DISYNI,
                      'G_isyn_dend': DGISYN 
                      }  

        # numerical integration
        k = 1
        start_time = time.time()    

        while r1.successful() and k < num_steps:
            r1.integrate(round(r1.t, 9) + self.t_dt) 
            # r1.integrate(r1.t + self.t_dt)
            T[k]=r1.t
            # if k < 10:
            #     print("num_steps ;", k, ", r1.t ;", r1.t) # CHJ-check  
            VS[k]=r1.y[0]
            SCA[k]=r1.y[1]
            SNAM[k]=r1.y[2]
            SNAH[k]=r1.y[3]
            SNAPM[k]=r1.y[4]
            SKDR[k]=r1.y[5]
            SCAM[k]=r1.y[6]
            SCAH[k]=r1.y[7]
            SHM[k]=r1.y[8]         
            VD[k]=r1.y[9]
            DCA[k]=r1.y[10]
            DCAL[k]=r1.y[11]
            DNAM[k]=r1.y[12]
            DNAH[k]=r1.y[13]
            DNAPM[k]=r1.y[14]
            DKDR[k]=r1.y[15]
            DKCAM[k]=r1.y[16]
            DCAM[k]=r1.y[17]
            DCAH[k]=r1.y[18]
            DHM[k]=r1.y[19]
            
            # Isoma
            if(self.is_sigType=='Step'):
                i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.is_heav_param
                IS[k] = i0 + s*((self.heav(poff1-r1.t)*self.heav(r1.t-pon1)*ip1)
                             + (self.heav(poff2-r1.t)*self.heav(r1.t-pon2)*ip2)
                             + (self.heav(poff3-r1.t)*self.heav(r1.t-pon3)*ip3)
                             + (self.heav(poff4-r1.t)*self.heav(r1.t-pon4)*ip4)
                             + (self.heav(poff5-r1.t)*self.heav(r1.t-pon5)*ip5))
            elif(self.is_sigType=='Ramp'):
                IS[k]=self.is_amp-((self.is_amp-self.is_iv)/self.is_period)*abs(r1.t-self.is_period)+self.is_offset
            elif(self.is_sigType=='Import'):
                IS[k]=np.interp(r1.t, self.is_time,  self.Is, 0, 0)
            elif(self.is_sigType=='sine'):
                IS[k]=self.is_amp*sin(2*pi*r1.t/self.is_period)+self.is_offset
            elif(self.is_sigType=='square'):
                x_t=r1.t%self.is_period
                if(x_t<(self.is_period/2)):
                    IS[k]=1*self.is_amp+self.is_offset
                elif(x_t>=self.is_period/2):
                    IS[k]=0+self.is_offset
            
            # [soma]
            # calcium reversal potential  
            if(SCA[k] < 1.e-100):
                SCA[k] = 1.e-100
            if(sEca != 0.):
                SVCA[k] = sEca
            else:
                SVCA[k]=(const*log(sCAo/SCA[k]))-70 # Vrest 보정
            # I_Can
            SCAI[k] = sgca*SCAM[k]**2*SCAH[k]*(VS[k]-SVCA[k])
            # I_Naf
            SNAI[k] = sgna*SNAM[k]**3*SNAH[k]*(VS[k]-svna)
            # I_Nap
            SNAPI[k] = sgnap*SNAPM[k]**3*(VS[k]-svna)
            # I_Kdr
            SKDRI[k] = sgkdr*SKDR[k]**4*(VS[k]-svk)
            # I_Kca
            SKCAI[k] = sgkca*(SCA[k]/(SCA[k]+skd))*(VS[k]-svk)
            # I_H
            SHI[k] = sgh*SHM[k]*(VS[k]-svh)
            # I_eSyn
            SGESYN[k] = np.interp(r1.t, self.sesyn_t, self.Se_syncon, 0, 0)
            SESYNI[k] = SGESYN[k]*(VS[k]-svesyn)
            # I_iSyn            
            SGISYN[k] = np.interp(r1.t, self.sisyn_t, self.Si_syncon, 0, 0)
            SISYNI[k] = SGISYN[k]*(VS[k]-svisyn)

            # [dendrite]
            # calcium reversal potential
            if(DCA[k] < 1.e-100):
                DCA[k] = 1.e-100
            if(dEca != 0.):
                DVCA[k] = dEca
            else:
                DVCA[k]=(const*log(dCAo/DCA[k]))-70 # Vrest 보정
            # I_Cal
            DCALI[k] = dgcal*DCAL[k]*(VD[k]-DVCA[k])
            # I_Naf
            DNAI[k] = dgna*DNAM[k]**3*DNAH[k]*(VD[k]-dvna)
            # I_Nap
            DNAPI[k] = dgnap*DNAPM[k]**3*(VD[k]-dvna)
            # I_Kdr
            DKDRI[k] = dgkdr*DKDR[k]**4*(VD[k]-dvk)
            # I_Kca
            DKCAI[k] = dgkca*DKCAM[k]*(VD[k]-dvk)
            # I_Can
            DCAI[k] = dgca*DCAM[k]**2*DCAH[k]*(VD[k]-DVCA[k])
            # I_H
            DHI[k] = dgh*DHM[k]*(VD[k]-dvh)
            # I_eSyn
            DGESYN[k] = np.interp(r1.t, self.desyn_t, self.De_syncon, 0, 0)
            DESYNI[k] = DGESYN[k]*(VD[k]-dvesyn)
            # I_iSyn
            DGISYN[k] = np.interp(r1.t, self.disyn_t, self.Di_syncon, 0, 0)
            DISYNI[k] = DGISYN[k]*(VD[k]-dvisyn)
            
            # spike detection & firing rate
            SpikeTimes, S_Detect = detect_Spike(SpikeTimes, r1.t, VS[k])
            if(S_Detect == True):
                cal_FR = cal_FiringRate(SpikeTimes)
                dt = self.t_dt
                i = int(round(r1.t/dt, 9))
                FR[i] = cal_FR
                FiringRate.append(cal_FR)
            
            k += 1

        # realtime taken for simulation     
        simulTime=time.time() - start_time
        # print "Elapsed time: ", simulTime, "s" # CHJ
        print("Elapsed time: ", simulTime, "s")
        
        self.resultArray=resultArray
        return resultArray
        
    # intracellular current injection
    def setIsSignal(self, signalType, period, amplitude, ivalue, offset, heav_param, time, Is):        
        self.is_sigType = signalType
        self.is_period = period
        self.is_amp = amplitude
        self.is_iv = ivalue
        self.is_offset = offset
        self.is_heav_param = heav_param
        self.is_time = time
        self.Is = Is
        
        # print "Intracellular current input set." # CHJ
        print("Intracellular current input set.")
    
    # synaptic conductance variation
    def setSynConSignal(self, G_e, G_i, e_times, i_times, department):
        if(department=='Soma'):
            self.Se_syncon=G_e
            self.Si_syncon=G_i
            self.sesyn_t=e_times
            self.sisyn_t=i_times
        elif(department=='Dendrite'):
            self.De_syncon=G_e
            self.Di_syncon=G_i
            self.desyn_t=e_times
            self.disyn_t=i_times
        
        # print "Synaptic conductance input set." # CHJ
        print("Synaptic conductance input set.")
        

## Child class of Cell class for muscle-tendon unit
class Musclefibers(Cell, object):
    def __init__(self, *args, **kwargs):
        super(Musclefibers, self).__init__(*args, **kwargs)
        self.setCellType('musclefibers')
        # print "Muscle-tendon unit instance constructed. ID: "+str(self.uniqueNumber) # CHJ
        print("Muscle-tendon unit instance constructed. ID: "+str(self.uniqueNumber))
        
    def __del__(self):
        # print "Muscle-tendon unit instance destructed. ID: "+str(self.uniqueNumber) # CHJ
        print("Muscle-tendon unit instance destructed. ID: "+str(self.uniqueNumber))
        
    def model(self, t, y):
        # state variables
        CaSR, CaSRCS, CaSP, CaSPB, CaSPT, C1, C2, A, XCE = y

        # model parameters
        K1, K2, K3, K4, K5i, K6i, K, Rmax, Umax, tau1, tau2,  phi1, phi2, phi3, phi4, CS0, B0, T0,  c1i,c1n1,c1n2,c1n3,c1n4, c2i,c2n1,c2n2,c2n3,c2n4, C3, C4, C5, alpha_i, a1, a2, a3, beta, gamma, KSE,  P0, g1, g2, g3, a0, b0, c0, d0, AM, LM, sf, cv = self.parameter
        ms=0.001
        b0, d0= b0*ms, d0*ms # unit conversion from s to ms
        
        # model equations
        n = len(y)      
        # dydt=range(n) # CHJ
        dydt = np.zeros(n)

        # [Module 1]
        # calcium release from sarcoplasmic reticulum
        R=self.cal_R(t, CaSR)
        
        # calcium uptake to sarcoplasmic reticulum
        U=Umax*(((CaSP**2)*(K**2))/(1+CaSP*K+(CaSP**2)*(K**2)))**2
     
        # length dependent activation under partial exciation
        Xm, Vm, Am=self.cal_Xm_Vm_Am(t)
        if(Xm<-8):
            uXm=(phi1*Xm)+phi2
        if(Xm>=-8):
            uXm=(phi3*Xm)+phi4
        K5=uXm*K5i
        
        # activation dependent calcium-troponin coupling
        K6=(K6i/(1+5*A))
     
        # d(CaSR)/dt
        dydt[0] = -K1*CS0*CaSR+(K1*CaSR+K2)*CaSRCS-R+U
        
        # CaSRCS
        dydt[1] = K1*CS0*CaSR-(K1*CaSR+K2)*CaSRCS
        
        # d(CaSP)/dt
        dydt[2] = -K5*T0*CaSP+(K5*CaSP+K6)*CaSPT-K3*B0*CaSP+(K3*CaSP+K4)*CaSPB+R-U
                 
        # d(CaSPB)/dt
        dydt[3] = K3*B0*CaSP-(K3*CaSP+K4)*CaSPB
        
        # d(CaSPT)/dt
        dydt[4] = K5*T0*CaSP-(K5*CaSP+K6)*CaSPT
        
        # [Module 2]
        # d(C1)/dt
        c1inf = c1n1*(1+tanh(((CaSPT/T0)-c1n2)/c1n3))+c1i
        c1tau = c1n4
        dydt[5] = (c1inf-C1)/c1tau
        
        # d(C2)/dt
        c2inf = c2n1*(1+tanh(((CaSPT/T0)-c2n2)/c2n3))+c2i
        c2tau = c2n4
        dydt[6] = (c2inf-C2)/c2tau
        
        # d(A_tilde)/dt
        Ainf=0.5*(1+tanh(((CaSPT/T0)-C1)/C2))
        TA=C3/(cosh(((CaSPT/T0)-C4)/(2*C5)))
        dydt[7]=(Ainf-A)/TA
        
        # A(t)
        As=self.cal_As(Xm,A,Vm,Am,t,alpha_i,a1,a2,a3,beta,gamma)
        
        # [Module 3]
        # length-tension property
        gXm=exp(-((Xm-g1)/g2)**2)+g3
        
        # velocity-tension property
        Fc=P0*gXm*As
        
        # force calculation
        F=self.cal_F(t, XCE, Xm)     
        
        # d(XCE)/dt
        if(F<=Fc):
            dydt[8]=(-b0*(Fc-F))/(F+P0*a0*Fc/P0)
        else:
            gainlength=(-d0*(Fc-F))/(2*Fc-F+P0*c0*Fc/P0)
            if(gainlength<=0):
                dydt[8]=(1.e-3)*(1.e5)
            else:
                dydt[8]=gainlength
  
        return dydt

    def solModel(self):
        # model parameters
        K1, K2, K3, K4, K5i, K6i, K, Rmax, Umax, tau1, tau2,  phi1, phi2, phi3, phi4, CS0, B0, T0,  c1i,c1n1,c1n2,c1n3,c1n4, c2i,c2n1,c2n2,c2n3,c2n4, C3, C4, C5, alpha_i, a1, a2, a3, beta, gamma, KSE,  P0, g1, g2, g3, a0, b0, c0, d0, AM, LM, sf, cv = self.parameter
        
        # numerical solver setting
        r2 = integrate.ode(self.model).set_integrator('lsoda', max_step=10*self.t_dt)
        r2.set_initial_value(self.ivalue, self.t_start)
        
        # arrays for integration results
        num_steps = int(np.floor((self.t_stop - self.t_start)/self.t_dt) + 1)  
        T2= np.zeros(num_steps)
        R=np.zeros(num_steps)
        F = np.zeros(num_steps)
        CaSR=np.zeros(num_steps)
        CaSRCS=np.zeros(num_steps)
        CaSP=np.zeros(num_steps)
        CaSPB=np.zeros(num_steps)
        CaSPT=np.zeros(num_steps)
        C1=np.zeros(num_steps)
        C2=np.zeros(num_steps)
        A=np.zeros(num_steps)
        XCE=np.zeros(num_steps)
        XM=np.zeros(num_steps)
        VM=np.zeros(num_steps)
        AM=np.zeros(num_steps)
        AS=np.zeros(num_steps)
        SP=np.zeros(num_steps)
        MUAP=np.zeros(num_steps) 

        # initialization of state variables 
        for i in self.sp_signal:
            k=int(round(i/(self.t_dt)))
            SP[k]=1
                
        CaSR[0]=self.ivalue[0]
        CaSRCS[0]=self.ivalue[1]
        CaSP[0]=self.ivalue[2]
        CaSPB[0]=self.ivalue[3]
        CaSPT[0]=self.ivalue[4]
        C1[0]=self.ivalue[5]
        C2[0]=self.ivalue[6]
        A[0]=self.ivalue[7]
        XCE[0]=self.ivalue[8]
        
        XM[0]=self.xm[0]
        F[0]=self.cal_F(T2[0], XCE[0], -8)
        AS[0]=self.cal_As(XM[0],A[0],VM[0], AM[0],T2[0],alpha_i,a1,a2,a3,beta,gamma)
        
        # numerical integration
        k = 1
        start_time = time.time()
        
        while r2.successful() and k < num_steps:
            r2.integrate(round(r2.t, 9) + self.t_dt)
            T2[k]=r2.t
            CaSR[k]=r2.y[0]
            CaSRCS[k]=r2.y[1]
            CaSP[k]=r2.y[2]
            CaSPB[k]=r2.y[3]
            CaSPT[k]=r2.y[4]
            C1[k]=r2.y[5]
            C2[k]=r2.y[6]
            A[k]=r2.y[7]
            XCE[k]=r2.y[8]
            
            R[k]=self.cal_R(T2[k], CaSR[k])
            XM[k], VM[k], AM[k]=self.cal_Xm_Vm_Am(T2[k])
            AS[k]=self.cal_As(XM[k],A[k],VM[k],AM[k],T2[k], alpha_i,a1,a2,a3,beta,gamma)
            MUAP[k]=self.cal_MUAP(T2[k])
            F[k]=self.cal_F(T2[k], XCE[k], XM[k]) 

            k += 1

        # realtime taken for simulation  
        self.simulTime=time.time() - start_time
        # print "Elapsed time: ", time.time() - start_time, "s" # CHJ  
        print("Elapsed time: ", time.time() - start_time, "s")
        
        # array for variables that are saved in result files
        resultArray2={ 'Time': T2,
                        'Xm': XM,
                       'Vm'  : VM,
                       'Am' : AM,
                       'R'  : R,
                       'Spike': SP,
                        'F': F,
                        'A_tilde': A,
                      'CaSR': CaSR,
                      'CaSRCS': CaSRCS,
                      'CaSP': CaSP,
                      'CaSPB': CaSPB,
                      'CaSPT': CaSPT,
                      'C1': C1,
                      'C2' : C2,
                      'A': AS,
                      'XCE': XCE,
                      'MUAP': MUAP
                      }

        self.resultArray= resultArray2
        
        return resultArray2      

    # Muscle-tendon unit length variation
    def setXmSignal(self, signalType, times, xm):

        self.XmSignalType = signalType
        self.xm_index = times
        self.xm = xm
        Xm = self.xm
        lengh = len(Xm)
        dt = self.t_dt 
        ms=0.001
        
        self.vm = np.zeros(lengh)
        self.am = np.zeros(lengh)
        self.vm[0] = 0.
        self.am[0] = 0.

        # calculation of velocity and acceleration using euler method
        if(self.XmSignalType=='Isometric'):
            pass

        else:
            if(self.XmSignalType=='Exp' or self.XmSignalType=='Random'):
                dt = self.xm_index[1] - self.xm_index[0]
                    
            # for t in xrange(1, lengh): # CHJ
            for t in range(1, lengh):
                xm = Xm[t]

                if(t == 1):
                    xm_1 = Xm[0]
                    self.vm[t] = (xm-xm_1)/(dt*ms)                
                    self.am[t] = self.vm[t]/dt/ms
                    self.am[t] = 0
                    
                else:
                    xm_1 = Xm[t-1]
                    xm_2 = Xm[t-2]
                    self.vm[t] = (xm-xm_1)/(dt*ms)
                    vm_1 = (xm_1-xm_2)/(dt*ms)
                    self.am[t] = (self.vm[t]-vm_1)/(dt*ms)
                
                if(abs(self.am[t]) < 1.e-3):
                    self.am[t] = 0.
                    
        # print "Muscle-tendon unit length set." # CHJ
        print("Muscle-tendon unit length set.")

    # current impulse stimulation with axonal conduction delay
    def setSpikeSignal(self, array):
        cv = self.parameter[49]
        ms = 1000
        sd = np.round(1/cv*ms, 1)
        self.sp_signal=array + sd
        
        # fitting to simulation time
        bool_arr = self.sp_signal <= self.t_stop
        self.sp_signal = self.sp_signal[bool_arr]

    # length, velocity, and acceleration
    def cal_Xm_Vm_Am(self, t):
       # interpolating Xm, Vm, Am
        if(self.XmSignalType=='Isometric'):
            xm=self.xm[0]
            vm=0.
            am=0.
        elif(self.XmSignalType=='Isokinetic'):
            xm=np.interp(t, self.xm_index,  self.xm, -8, self.xm[-1])
            vm=np.interp(t, self.xm_index,  self.vm, -1, -8)
            am=np.interp(t, self.xm_index,  self.am, -1, -8)
                
        elif(self.XmSignalType=='Random'):
            xm=np.interp(t, self.xm_index,  self.xm)
            vm=np.interp(t, self.xm_index,  self.vm)
            am=np.interp(t, self.xm_index,  self.am)
                
        elif(self.XmSignalType=='Exp'):
            xm=np.interp(t, self.xm_index,  self.xm, -8, self.xm[-1])
            vm=np.interp(t, self.xm_index,  self.vm, -1, -8)
            am=np.interp(t, self.xm_index,  self.am, -1, -8)
     
        return xm, vm, am   
    
    # calcium release from sarcoplasmic reticulum
    def cal_R(self, t ,CaSR):
        R=0.
        Rmax=self.parameter[7]
        tau1, tau2=self.parameter[9], self.parameter[10]
 
        for ts in self.sp_signal:
            if(t>=ts):
                # 2020.03.13 revision : Pmax -> Rmax
                R+=CaSR*Rmax*(1-exp(-(t-ts)/tau1))*exp(-(t-ts)/tau2)
        
        return R

    # muscle unit action potential
    def cal_MUAP(self, t):
        AM = self.parameter[46]
        LM = self.parameter[47]
        HR1 = 0.
        HR2 = 0.

        for ts in self.sp_signal:
            HR1+=AM*(t-ts)*np.exp(-((t-ts)/LM)**2)
            HR2+=AM*(1-2*(((t-ts)/LM)**2))*np.exp(-((t-ts)/LM)**2)

        return HR1/2 + HR2

    # muscle activation level
    def cal_As(self, Xm,A,Vm,Am,t, alpha_i,a1,a2,a3,beta,gamma):
        L0=-8
        phi1, phi2, phi3, phi4 = self.parameter[11:15]

        if(Xm<L0):
            uxm=(phi1*Xm)+phi2
        if(Xm>=L0):
            uxm=(phi3*Xm)+phi4
        
        # detection of dynamic length change
        if(Am == 0.):
            As=A**alpha_i
        else:    
            alpha=a1*(1+tanh((t-a2)/a3))+alpha_i
            At=A**(alpha)
            if((Xm<=L0) and (Vm>=0)):
                Ad=(1+beta*uxm)*(1+gamma*Vm)
                As=At/Ad
            else:
                As=At

        if(As<1.e-3):
            As=1.e-3
           
        return As

    # force production
    def cal_F(self, t, XCE, Xm):
        xm_init=-8
        xce_init=-8
        
        d_xce=XCE-xce_init
        d_xm=Xm-xm_init
        d_se=d_xm-d_xce
        
        P0=self.parameter[38]
        KSE=self.parameter[37]
        sf=self.parameter[48]

        if(t<0.002): 
            F=1.e-5
        elif(d_se<=0):
            F=0.
        else:
            F=P0*KSE*(d_se)
            
        F=sf*F

        return F
    

## Parent class for a unit
class Unit(object):
    def __init__(self, *args, **kwargs):
        self.setID(kwargs['uniqueNumber'])
        # print "Unit instance constructed. ID: "+str(self.uniqueNumber) # CHJ
        print("Unit instance constructed. ID: "+str(self.uniqueNumber))
        
    def __del__(self):
        # print "Unit instance destructed. ID: "+str(self.uniqueNumber) # CHJ
        print("Unit instance destructed. ID: "+str(self.uniqueNumber))

    def setUnitType(self, unitType):
        self.unitType=unitType

    def getUnitType(self):
        return self.unitType

    def setID(self, num):
        self.uniqueNumber=num

    def getID(self):
        return self.uniqueNumber
        
    def getHostname(self):
        self.hostname = platform.node()
        return self.hostname

    def setSimulTime(self, t_start, t_stop, resolution):
        self.t_start=t_start
        self.t_stop=t_stop
        self.t_dt=resolution

        # print "Simulation time changed. ID: "+str(self.getID()) # CHJ
        print("Simulation time changed. ID: "+str(self.getID()))
        # print "start time: "+str(self.t_start)+" (msec)" # CHJ
        print("start time: "+str(self.t_start)+" (msec)")
        # print "stop time: "+str(self.t_stop)+" (msec)" # CHJ
        print("stop time: "+str(self.t_stop)+" (msec)")
        # print "time interval: "+str(self.t_dt)+" (msec)" # CHJ
        print("time interval: "+str(self.t_dt)+" (msec)")

    def getSimulTime(self):
        return self.t_stop, self.t_dt
        
    def setInitialValue(self, ivalue_1, ivalue_2):
        self.m_ivalue=ivalue_1
        self.f_ivalue=ivalue_2
        # print str(self.getUnitType())+" initial values set." # CHJ
        print(str(self.getUnitType())+" initial values set.")

    def getInitialValue(self):
        return self.ivalue


## Child class of Unit class for motor unit
class Motorunit(Unit, object):
    def __init__(self, *args, **kwargs):
        super(Motorunit, self).__init__(*args, **kwargs)
        self.setUnitType('motorunit')
        self.uniqueNumber=kwargs['uniqueNumber']
        self.MN = Motoneuron(uniqueNumber=kwargs['uniqueNumber'])
        self.MC = Musclefibers(uniqueNumber=kwargs['uniqueNumber'])
        
        # print "Motorunit instance constructed. ID: "+str(kwargs['uniqueNumber']) # CHJ
        print("Motorunit instance constructed. ID: "+str(kwargs['uniqueNumber']))

    def __del__(self):
        # print "Motorunit instance destructed. ID: "+str(self.uniqueNumber) # CHJ
        print("Motorunit instance destructed. ID: "+str(self.uniqueNumber))

    def setParameter(self, param_1, param_2):
        self.m_parameter=param_1
        self.f_parameter=param_2
        
        self.MN.setParameter(param_1)
        self.MC.setParameter(param_2)

    def getParameter(self):
        return self.m_parameter, self.f_parameter

    def setSimulTime(self, t_start, t_stop, resolution):
        Unit.setSimulTime(self, t_start, t_stop, resolution)
        self.MN.setSimulTime(t_start, t_stop, resolution)
        self.MC.setSimulTime(t_start, t_stop, resolution)

    def setInitialValue(self, ivalue_1, ivalue_2):
        Unit.setInitialValue(self, ivalue_1, ivalue_2)
        self.MN.setInitialValue(ivalue_1)
        self.MC.setInitialValue(ivalue_2)

    def setIsSignal(self, signalType, period, amplitude, ivalue, offset, heav_param, time, Is):
        self.MN.setIsSignal(signalType, period, amplitude, ivalue, offset, heav_param, time, Is)
       
        # print "Intracellular current input set." # CHJ
        print("Intracellular current input set.")
        
    def setSynConSignal(self, G_e, G_i, e_times, i_times, department):
        self.MN.setSynConSignal(G_e, G_i, e_times, i_times, department)
                
        # print "Synaptic conductance input set." # CHJ
        print("Synaptic conductance input set.")
        
    def setXmSignal(self, signalType, times, xm):
        self.MC.setXmSignal(signalType, times, xm)

    def solModel(self):
        # model parameters for motoneuron
        gms,gmd,gc,cmd,cms, \
        svl,dvl, \
        sf,skca,salpha,sCAo,sEca, \
        sgna,svna,sanamc,sanamv,sanama,sanamb,sbnamc,sbnamv,sbnama,sbnamb,snahth,snahslp,snahv,snaha,snahb,snahc, \
        sgnap,svna,sanapmc,sanapmv,sanapma,sanapmb,sbnapmc,sbnapmv,sbnapma,sbnapmb, \
        sgkdr,svk,skdrth,skdrslp,skdrv,skdra,skdrb,skdrc, \
        sgkca,svk,skd, \
        sgca,scamth,scamslp,scamtau,scahth,scahslp,scahtau, \
        sgh,svh,shth,shslp,shtau, \
        svesyn, \
        svisyn, \
        df,dkca,dalpha,dCAo,dEca, \
        dgcal,dcalth,dcalslp,dcaltau, \
        dgna,dvna,danamc,danamv,danama,danamb,dbnamc,dbnamv,dbnama,dbnamb,dnahth,dnahslp,dnahv,dnaha,dnahb,dnahc, \
        dgnap,dvna,danapmc,danapmv,danapma,danapmb,dbnapmc,dbnapmv,dbnapma,dbnapmb, \
        dgkdr,dvk,dkdrth,dkdrslp,dkdrv,dkdra,dkdrb,dkdrc, \
        dgkca,dvk,dkcamth,dkcamslp,dkcamtau, \
        dgca,dcamth,dcamslp,dcamtau,dcahth,dcahslp,dcahtau, \
        dgh,dvh,dhth,dhslp,dhtau, \
        dvesyn, \
        dvisyn=self.MN.parameter
        
        # model parameters for muscle-tendon unit
        K1, K2, K3, K4, K5i, K6i, K, Rmax, Umax, tau1, tau2,  phi1, phi2, phi3, phi4, CS0, B0, T0,  c1i,c1n1,c1n2,c1n3,c1n4, c2i,c2n1,c2n2,c2n3,c2n4, C3, C4, C5, alpha_i, a1, a2, a3, beta, gamma, KSE,  P0, g1, g2, g3, a0, b0, c0, d0, AM, LM, sf, cv = self.MC.parameter

        # constants for motoneuron
        R=8.31441
        Temp=309.15
        Zca=2
        Fe=96485.309
        const = 1000*R*Temp/Zca/Fe
        
        SpikeTimes=[]

        # setting of numerical integration for motoneuron
        r1 = integrate.ode(self.MN.model).set_integrator('vode')
        #r1 = integrate.ode(self.MN.model).set_integrator('lsoda', max_step=10*self.t_dt)
        r1.set_initial_value(self.MN.ivalue, self.t_start)
                
        num_steps = int(np.floor((self.t_stop - self.t_start)/self.t_dt) + 1)
        FR = np.zeros(num_steps)
        T = np.zeros(num_steps)
        IS = np.zeros(num_steps)
        VS = np.zeros(num_steps)
        SCA = np.zeros(num_steps)
        SVCA = np.zeros(num_steps)
        SNAI = np.zeros(num_steps)
        SNAM = np.zeros(num_steps)
        SNAH = np.zeros(num_steps)
        SNAPI = np.zeros(num_steps)
        SNAPM = np.zeros(num_steps) 
        SKDRI = np.zeros(num_steps)
        SKDR = np.zeros(num_steps)
        SKCAI = np.zeros(num_steps)
        SCAI = np.zeros(num_steps)
        SCAM = np.zeros(num_steps) 
        SCAH = np.zeros(num_steps)
        SHI = np.zeros(num_steps)
        SHM = np.zeros(num_steps)
        VD = np.zeros(num_steps) 
        DCA = np.zeros(num_steps)
        DVCA = np.zeros(num_steps)
        DCALI = np.zeros(num_steps)
        DCAL = np.zeros(num_steps) 
        DNAI = np.zeros(num_steps)
        DNAM = np.zeros(num_steps)
        DNAH = np.zeros(num_steps)
        DNAPI = np.zeros(num_steps)
        DNAPM = np.zeros(num_steps)
        DKDRI = np.zeros(num_steps)
        DKDR = np.zeros(num_steps)
        DKCAI = np.zeros(num_steps)
        DKCAM = np.zeros(num_steps)
        DCAI = np.zeros(num_steps)
        DCAM = np.zeros(num_steps)
        DCAH = np.zeros(num_steps)
        DHI = np.zeros(num_steps)
        DHM = np.zeros(num_steps)
        SESYNI = np.zeros(num_steps)
        SGESYN = np.zeros(num_steps)
        SISYNI = np.zeros(num_steps)
        SGISYN = np.zeros(num_steps)
        DESYNI = np.zeros(num_steps)
        DGESYN = np.zeros(num_steps)
        DISYNI = np.zeros(num_steps)
        DGISYN = np.zeros(num_steps)
        MN_SP = np.zeros(num_steps)

        IS[0]= self.MN.Is[0]
        VS[0]=self.MN.ivalue[0]
        SCA[0]=self.MN.ivalue[1]
        if(sEca != 0.):
            SVCA[0] = sEca
        else:
            SVCA[0]=(const*log(sCAo/SCA[0]))-70 # Vrest 보정
        SNAM[0]=self.MN.ivalue[2]
        SNAH[0]=self.MN.ivalue[3]
        SNAI[0] = sgna*SNAM[0]**3*SNAH[0]*(VS[0]-svna)
        SNAPM[0]=self.MN.ivalue[4]
        SNAPI[0] = sgnap*SNAPM[0]**3*(VS[0]-svna)
        SKDR[0]=self.MN.ivalue[5]
        SKDRI[0] = sgkdr*SKDR[0]**4*(VS[0]-svk)
        SKCAI[0] = sgkca*(SCA[0]/(SCA[0]+skd))*(VS[0]-svk)
        SCAM[0]=self.MN.ivalue[6]
        SCAH[0]=self.MN.ivalue[7]
        SCAI[0] = sgca*SCAM[0]**2*SCAH[0]*(VS[0]-SVCA[0])
        SHM[0]=self.MN.ivalue[8]
        SHI[0] = sgh*SHM[0]*(VS[0]-svh)
        VD[0]=self.MN.ivalue[9]
        DCA[0]=self.MN.ivalue[10]
        if(dEca != 0.):
            DVCA[0] = dEca
        else:
            DVCA[0]=(const*log(dCAo/DCA[0]))-70
        DCAL[0]=self.MN.ivalue[11]
        DCALI[0] = dgcal*DCAL[0]*(VD[0]-DVCA[0])
        DNAM[0]=self.MN.ivalue[12]
        DNAH[0]=self.MN.ivalue[13]
        DNAI[0] = dgna*DNAM[0]**3*DNAH[0]*(VD[0]-dvna)
        DNAPM[0]=self.MN.ivalue[14]
        DNAPI[0] = dgnap*DNAPM[0]**3*(VD[0]-dvna)
        DKDR[0]=self.MN.ivalue[15]
        DKDRI[0] = dgkdr*DKDR[0]**4*(VD[0]-dvk)        
        DKCAM[0] = self.MN.ivalue[16]
        DKCAI[0] = dgkca*DKCAM[0]*(VD[0]-dvk)
        DCAM[0]=self.MN.ivalue[17]
        DCAH[0]=self.MN.ivalue[18]
        DCAI[0] = dgca*DCAM[0]**2*DCAH[0]*(VD[0]-DVCA[0])
        DHM[0]=self.MN.ivalue[19]
        DHI[0] = dgh*DHM[0]*(VD[0]-dvh)
        SGESYN[0] = self.MN.Se_syncon[0]
        SESYNI[0] = SGESYN[0]*(VS[0]-svesyn)
        SGISYN[0] = self.MN.Si_syncon[0]
        SISYNI[0] = SGISYN[0]*(VS[0]-svisyn)
        DGESYN[0] = self.MN.De_syncon[0]
        DESYNI[0] = DGESYN[0]*(VD[0]-dvesyn)
        DGISYN[0] = self.MN.Di_syncon[0]
        DISYNI[0] = DGISYN[0]*(VD[0]-dvisyn)
        
        resultArray={'Time': T,
                      'Is' : IS,
                      'Firing_rate' : FR,
                      'MN_Spike': MN_SP,
                      'V_soma': VS,
                      '[Ca]_soma': SCA,
                      'E_Ca_soma': SVCA,
                      'I_Naf_soma': SNAI,
                      'm_Naf_soma': SNAM,
                      'h_Naf_soma': SNAH,
                      'I_Nap_soma': SNAPI,
                      'm_Nap_soma': SNAPM,
                      'I_Kdr_soma': SKDRI,
                      'n_Kdr_soma' : SKDR,
                      'I_Kca_soma': SKCAI,
                      'I_Can_soma': SCAI,
                      'm_Can_soma': SCAM,
                      'h_Can_soma': SCAH,
                      'I_H_soma': SHI,
                      'm_H_soma' : SHM,    
                      'V_dend': VD,
                      '[Ca]_dend': DCA,
                      'E_Ca_dend': DVCA,
                      'I_Cal_dend': DCALI,
                      'm_Cal_dend': DCAL,
                      'I_Naf_dend': DNAI,
                      'm_Naf_dend': DNAM,
                      'h_Naf_dend': DNAH,
                      'I_Nap_dend': DNAPI,
                      'm_Nap_dend': DNAPM,
                      'I_Kdr_dend': DKDRI,
                      'n_Kdr_dend' : DKDR,
                      'I_Kca_dend': DKCAI,
                      'm_Kca_dend': DKCAM,
                      'I_Can_dend': DCAI,
                      'm_Can_dend': DCAM,
                      'h_Can_dend': DCAH,
                      'I_H_dend': DHI,
                      'm_H_dend' : DHM,
                      'I_esyn_soma': SESYNI,
                      'G_esyn_soma': SGESYN,
                      'I_isyn_soma': SISYNI,
                      'G_isyn_soma': SGISYN,
                      'I_esyn_dend': DESYNI,
                      'G_esyn_dend': DGESYN, 
                      'I_isyn_dend': DISYNI,
                      'G_isyn_dend': DGISYN 
                      }  
        
        # setting of numerical integration for muscle-tendon unit
        r2 = integrate.ode(self.MC.model).set_integrator('lsoda', rtol=1e-6, atol=1e-12, max_step=10*self.t_dt)
        #r2 = integrate.ode(self.MC.model).set_integrator('vode')
        r2.set_initial_value(self.MC.ivalue, self.t_start)
            
        num_steps = int(np.floor((self.t_stop - self.t_start)/self.t_dt) + 1)  
        T2= np.zeros(num_steps)
        R=np.zeros(num_steps)
        F = np.zeros(num_steps)
        CaSR=np.zeros(num_steps)
        CaSRCS=np.zeros(num_steps)
        CaSP=np.zeros(num_steps)
        CaSPB=np.zeros(num_steps)
        CaSPT=np.zeros(num_steps)
        C1=np.zeros(num_steps)
        C2=np.zeros(num_steps)
        A=np.zeros(num_steps)
        XCE=np.zeros(num_steps)
        XM=np.zeros(num_steps)
        VM=np.zeros(num_steps)
        AM=np.zeros(num_steps)
        AS=np.zeros(num_steps)
        MF_SP=np.zeros(num_steps)
        MUAP=np.zeros(num_steps) 

        CaSR[0]=self.MC.ivalue[0]
        CaSRCS[0]=self.MC.ivalue[1]
        CaSP[0]=self.MC.ivalue[2]
        CaSPB[0]=self.MC.ivalue[3]
        CaSPT[0]=self.MC.ivalue[4]
        C1[0]=self.MC.ivalue[5]
        C2[0]=self.MC.ivalue[6]
        A[0]=self.MC.ivalue[7]
        XCE[0]=self.MC.ivalue[8]
        XM[0]=self.MC.xm[0]
        F[0]=self.MC.cal_F(T2[0], XCE[0], -8)
        AS[0]=self.MC.cal_As(XM[0],A[0],VM[0], AM[0],T2[0],alpha_i,a1,a2,a3,beta,gamma)
        
        resultArray2={'Time': T2,
                      'Xm': XM,
                      'Vm'  : VM,
                      'Am' : AM,
                      'R'  : R,
                      'MF_Spike': MF_SP,
                      'F': F,
                      'A_tilde': A,
                      'CaSR': CaSR,
                      'CaSRCS': CaSRCS,
                      'CaSP': CaSP,
                      'CaSPB': CaSPB,
                      'CaSPT': CaSPT,
                      'C1': C1,
                      'C2' : C2,
                      'A': AS,
                      'XCE': XCE,
                      'MUAP': MUAP
                      }
                      
        resultArray.update(resultArray2)
        
        # numerical integration
        k = 1
        start_time = time.time()

        while r2.successful() and k < num_steps:
            # numerical integration for motoneuron
            r1.integrate(round(r1.t, 9) + self.t_dt)
            VS[k]=r1.y[0]
            SCA[k]=r1.y[1]
            SNAM[k]=r1.y[2]
            SNAH[k]=r1.y[3]
            SNAPM[k]=r1.y[4]
            SKDR[k]=r1.y[5]
            SCAM[k]=r1.y[6]
            SCAH[k]=r1.y[7]
            SHM[k]=r1.y[8]         
            VD[k]=r1.y[9]
            DCA[k]=r1.y[10]
            DCAL[k]=r1.y[11]
            DNAM[k]=r1.y[12]
            DNAH[k]=r1.y[13]
            DNAPM[k]=r1.y[14]
            DKDR[k]=r1.y[15]
            DKCAM[k]=r1.y[16]
            DCAM[k]=r1.y[17]
            DCAH[k]=r1.y[18]
            DHM[k]=r1.y[19]
            
            # Isoma
            if(self.MN.is_sigType=='Step'):
                i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.MN.is_heav_param
                IS[k] = i0 + s*((self.MN.heav(poff1-r1.t)*self.MN.heav(r1.t-pon1)*ip1)
                             + (self.MN.heav(poff2-r1.t)*self.MN.heav(r1.t-pon2)*ip2)
                             + (self.MN.heav(poff3-r1.t)*self.MN.heav(r1.t-pon3)*ip3)
                             + (self.MN.heav(poff4-r1.t)*self.MN.heav(r1.t-pon4)*ip4)
                             + (self.MN.heav(poff5-r1.t)*self.MN.heav(r1.t-pon5)*ip5))
            elif(self.MN.is_sigType=='Ramp'):
                IS[k]=self.MN.is_amp-((self.MN.is_amp-self.MN.is_iv)/self.MN.is_period)*abs(r1.t-self.MN.is_period)+self.MN.is_offset
            elif(self.MN.is_sigType=='Import'):
                IS[k]=np.interp(r1.t, self.MN.is_time,  self.MN.Is, 0, 0)
            elif(self.MN.is_sigType=='sine'):
                IS[k]=self.MN.is_amp*sin(2*pi*r1.t/self.MN.is_period)+self.MN.is_offset
            elif(self.MN.is_sigType=='square'):
                x_t=r1.t%self.MN.is_period
                if(x_t<(self.MN.is_period/2)):
                    IS[k]=1*self.MN.is_amp+self.MN.is_offset
                elif(x_t>=self.MN.is_period/2):
                    IS[k]=0+self.MN.is_offset
                    
            # [soma]
            # calcium reversal potential  
            if(SCA[k] < 1.e-100):
                SCA[k] = 1.e-100
            if(sEca != 0.):
                SVCA[k] = sEca
            else:
                SVCA[k]=(const*log(sCAo/SCA[k]))-70 # Vrest 보정
            # I_Can
            SCAI[k] = sgca*SCAM[k]**2*SCAH[k]*(VS[k]-SVCA[k])
            # I_Naf
            SNAI[k] = sgna*SNAM[k]**3*SNAH[k]*(VS[k]-svna)
            # I_Nap
            SNAPI[k] = sgnap*SNAPM[k]**3*(VS[k]-svna)
            # I_Kdr
            SKDRI[k] = sgkdr*SKDR[k]**4*(VS[k]-svk)
            # I_Kca
            SKCAI[k] = sgkca*(SCA[k]/(SCA[k]+skd))*(VS[k]-svk)
            # I_H
            SHI[k] = sgh*SHM[k]*(VS[k]-svh)
            # I_eSyn
            SGESYN[k] = np.interp(r1.t, self.MN.sesyn_t, self.MN.Se_syncon, 0, 0)
            SESYNI[k] = SGESYN[k]*(VS[k]-svesyn)
            # I_iSyn            
            SGISYN[k] = np.interp(r1.t, self.MN.sisyn_t, self.MN.Si_syncon, 0, 0)
            SISYNI[k] = SGISYN[k]*(VS[k]-svisyn)

            # [dendrite]
            # calcium reversal potential
            if(DCA[k] < 1.e-100):
                DCA[k] = 1.e-100
            if(dEca != 0.):
                DVCA[k] = dEca
            else:
                DVCA[k]=(const*log(dCAo/DCA[k]))-70
            # I_Cal
            DCALI[k] = dgcal*DCAL[k]*(VD[k]-DVCA[k])
            # I_Naf
            DNAI[k] = dgna*DNAM[k]**3*DNAH[k]*(VD[k]-dvna)
            # I_Nap
            DNAPI[k] = dgnap*DNAPM[k]**3*(VD[k]-dvna)
            # I_Kdr
            DKDRI[k] = dgkdr*DKDR[k]**4*(VD[k]-dvk)
            # I_Kca
            DKCAI[k] = dgkca*DKCAM[k]*(VD[k]-dvk)
            # I_Can
            DCAI[k] = dgca*DCAM[k]**2*DCAH[k]*(VD[k]-DVCA[k])
            # I_H
            DHI[k] = dgh*DHM[k]*(VD[k]-dvh)
            # I_eSyn
            DGESYN[k] = np.interp(r1.t, self.MN.desyn_t, self.MN.De_syncon, 0, 0)
            DESYNI[k] = DGESYN[k]*(VD[k]-dvesyn)
            # I_iSyn
            DGISYN[k] = np.interp(r1.t, self.MN.disyn_t, self.MN.Di_syncon, 0, 0)
            DISYNI[k] = DGISYN[k]*(VD[k]-dvisyn)
            # spike detection
            SpikeTimes, S_Detect = detect_Spike(SpikeTimes, r1.t, VS[k])

            if(S_Detect == True):
                cal_FR = cal_FiringRate(SpikeTimes)
                dt = self.t_dt
                i = int(round(r1.t/dt, 9))
                FR[i] = cal_FR
                MN_SP[i] = 1

            self.MC.setSpikeSignal(SpikeTimes)

            # numerical integration for muscle-tendon unit
            r2.integrate(round(r2.t, 9) + self.t_dt)
            T2[k]=r2.t
            CaSR[k]=r2.y[0]
            CaSRCS[k]=r2.y[1]
            CaSP[k]=r2.y[2]
            CaSPB[k]=r2.y[3]
            CaSPT[k]=r2.y[4]
            C1[k]=r2.y[5]
            C2[k]=r2.y[6]
            A[k]=r2.y[7]
            XCE[k]=r2.y[8]
            R[k]=self.MC.cal_R(T2[k], CaSR[k])
            XM[k], VM[k], AM[k]=self.MC.cal_Xm_Vm_Am(T2[k])
            AS[k]=self.MC.cal_As(XM[k],A[k],VM[k],AM[k],T2[k], alpha_i,a1,a2,a3,beta,gamma)
            MUAP[k]=self.MC.cal_MUAP(T2[k])
            F[k]=self.MC.cal_F(T2[k], XCE[k], XM[k])         
  
            k += 1

        # realtime taken for simulation
        self.simulTime=time.time() - start_time
        # print "Elapsed Time: ", time.time() - start_time, "s" # CHJ 
        print("Elapsed Time: ", time.time() - start_time, "s")
        
        # times of MUAP
        for i in self.MC.sp_signal:
            k=int(round(i/(self.t_dt)))
            MF_SP[k]=1

        self.resultArray=resultArray

        return resultArray

## Parent class of a population
class Pool(object):
    def __init__(self, *args, **kwargs):
        self.setUnitCount(kwargs['unit_count'])
        # print "Population instance constructed." # CHJ
        print("Population instance constructed.")

    def setUnitType(self, unitType):
        self.unit_type=unitType

    def getUnitType(self):
        return self.unit_type

    def setUnitsArray(self, array):
        self.arr_units=array

    def getUnitsArray(self):
        return self.arr_units

    def createCells(self, unit_count):
        self.setUnitsArray([]) 
        self.setUnitCount(unit_count)

    def setUnitCount(self, unit_count):
        self.unit_count=unit_count

    def getUnitCount(self):
        return self.unit_count

    def setSimulTimes(self, t_start, t_stop, t_dt):
        self.t_start=t_start
        self.t_stop=t_stop
        self.t_dt=t_dt
        arr_units=self.getUnitsArray()
        
        for unit in arr_units:
            unit.setSimulTime(t_start, t_stop, t_dt)

    def getSimulTimes(self):
        return self.t_start, self.t_stop, self.t_dt
        
    def setParameters(self, params):
        for i in range(self.unit_count):
            self.arr_units[i].setParameter(params[i])

    def setInitialValues(self, ivalue):
        arr_units=self.getUnitsArray()
        
        for unit in arr_units:
            unit.setInitialValue(ivalue)        
        
    def getInitialValues(self):
        pass


## Child class of Pool class for motoneuron population        
class MotoneuronPool(Pool, object):
    def __init__(self, *args, **kwargs):
        Pool.__init__(self, *args, **kwargs)
        self.setUnitType('motoneuron')
        self.createCells(self.unit_count)

    def createCells(self, unit_count):
        Pool.createCells(self, unit_count)
        arr_units=self.getUnitsArray()
        
        for i in range(1, unit_count+1):
            cell=Motoneuron(uniqueNumber=i)
            arr_units.append(cell)   
        
        self.setUnitsArray(arr_units)
   
    def setNeuronInputSignals(self, signalType, period, amplitude, ivalue, offset, heav_param, time, Is, MN_range_1, MN_range_2):
        arr_units=self.getUnitsArray()

        for i, unit in enumerate(arr_units):
            i = i+1
            if(i>=MN_range_1 and i <= MN_range_2):
                unit.setIsSignal(signalType, period, amplitude, ivalue, offset, heav_param, time.copy(), Is.copy())
            
    def setSynConSignals(self, G_e, G_i, e_times, i_times, MN_range_1, MN_range_2, department):
        arr_units=self.getUnitsArray() 
        
        for i, unit in enumerate(arr_units):
            i = i+1
            if(i>=MN_range_1 and i <= MN_range_2):
                unit.setSynConSignal(G_e.copy(), G_i.copy(), e_times.copy(), i_times.copy(), department)

    def getNeuronInputSignals(self):
        pass
    

## Child class of Pool class for muscle-tendon unit population    
class MusclefibersPool(Pool, object):
    def __init__(self, *args, **kwargs):
        Pool.__init__(self, *args, **kwargs)
        self.setUnitType('musclefibers')
        self.createCells(self.unit_count)

    def createCells(self, unit_count):
        Pool.createCells(self, unit_count)
        arr_units=self.getUnitsArray()
        
        for i in range(1, unit_count+1):
            cell=Musclefibers(uniqueNumber=i)
            arr_units.append(cell)   
        
        self.setUnitsArray(arr_units)

    def setSpikeSignals(self, sp_signal, MF_range_1, MF_range_2):
        arr_units=self.getUnitsArray() 
        
        for i, unit in enumerate(arr_units):
            i = i+1
            if(i>=MF_range_1 and i <= MF_range_2):
                unit.setSpikeSignal(sp_signal.copy())
                    
    def getSpikeSignals(self):
        pass

    def setMuscleLengthSignals(self, signalType, times, xm, MF_range_1, MF_range_2):
        arr_units=self.getUnitsArray() 
        
        for i, unit in enumerate(arr_units):
            i = i+1
            if(i>=MF_range_1 and i <= MF_range_2):
                unit.setXmSignal(signalType, times.copy(), xm.copy())        
                
    def getMuscleLengthSignals(self):
        pass


## Child class of Pool class for motor unit population    
class MotorunitPool(Pool, object):
    def __init__(self, *args, **kwargs):
        Pool.__init__(self, *args, **kwargs)
        self.setUnitType('motorunit')
        self.createCells(self.unit_count)

    def createCells(self, unit_count):
        Pool.createCells(self, unit_count)
        arr_units=self.getUnitsArray()
        
        for i in range(1, unit_count+1):
            unit=Motorunit(uniqueNumber=i)
            arr_units.append(unit)   
        
        self.setUnitsArray(arr_units)

    def setParameters(self, params_1, params_2):
        for i in range(self.unit_count):
            self.arr_units[i].setParameter(params_1[i], params_2[i])

    def setInitialValues(self, ivalue_1, ivalue_2):
        arr_units=self.getUnitsArray()
        
        for unit in arr_units:
            unit.setInitialValue(ivalue_1, ivalue_2)   

    def setNeuronInputSignals(self, signalType, period, amplitude, ivalue, offset, heav_param, time, Is, MN_range_1, MN_range_2):
        arr_units=self.getUnitsArray() 
        
        for i, unit in enumerate(arr_units):
            i = i+1
            if(i>=MN_range_1 and i <= MN_range_2):
                unit.setIsSignal(signalType, period, amplitude, ivalue, offset, heav_param, time.copy(), Is.copy())
        
    def getNeuronInputSignals(self):
        pass
    
    def setSynConSignals(self, G_e, G_i, e_times, i_times, MN_range_1, MN_range_2, department):
        arr_units=self.getUnitsArray() 
        
        for i, unit in enumerate(arr_units):
            i = i+1
            if(i>=MN_range_1 and i <= MN_range_2):
                unit.setSynConSignal(G_e.copy(), G_i.copy(), e_times.copy(), i_times.copy(), department)

    def setMuscleLengthSignals(self, signalType, times, xm, MF_range_1, MF_range_2):
        arr_units=self.getUnitsArray() 
        
        for i, unit in enumerate(arr_units):
            i = i+1
            if(i>=MF_range_1 and i <= MF_range_2):
                unit.setXmSignal(signalType, times.copy(), xm.copy())       
                
            
## Association class for intracellular current stimulation for motoneurons         
class NeuronSignalGenerator:
    def __init__(self, MP, unitCount, signalType='Ramp', period=5000., amplitude=0., offset=0., ivalue=0., time=10000., dt=0.1, heav_param=[]):
        self.genSignal(signalType, period, amplitude, offset, ivalue, time, dt, heav_param)
        self.setSignal(MP, 1, unitCount)
        
    def heav(self, x):
        return np.where(x < 0, 0, 1)
        
    def importSignalFile(self, filePath):
        self.signalType = 'Import'
        data=pd.read_csv(filePath)
        self.time = data.time
        self.Is = data.Is

    def genSignal(self, signalType, period, amplitude, offset, ivalue, time, dt, heav_param):
        p=float(period)
        a=float(amplitude)
        iv=float(ivalue)
        o=float(offset)
        
        num = int((time)/dt)+1
        t = np.linspace(0, time, num, endpoint=True)
        Is = np.zeros(num)
        k = 0
        
        if(signalType=='Step'):
            i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = heav_param
            Is = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                         + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                         + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                         + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                         + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
                         
        elif(signalType=='Ramp'):
            Is = a-((a-iv)/p)*np.abs(t-p)+o

        elif(signalType=='sine'):
            Is = a*np.sin(2*np.pi*t/p)+o #Is[k] -> Is
                
        elif(signalType=='square'):
            for i in t:
                x_t=i%period
                if(x_t<(period/2)):
                    Is[k]=1*a+o
                elif(x_t>=period/2):
                    Is[k]=0+o
                k+=1

        self.time = t
        self.Is = Is
        self.signalType = signalType
        self.period = period
        self.amplitude = amplitude
        self.ivalue = ivalue
        self.offset = offset
        self.heav_param = heav_param

    def setSignal(self, MP, MN_range_1, MN_range_2):
        MP.setNeuronInputSignals(self.signalType, self.period, self.amplitude, self.ivalue, self.offset, self.heav_param, self.time, self.Is, MN_range_1, MN_range_2)
        
    def plotSignal(self):
        fig = plt.figure(dpi=100)
        ax=fig.add_subplot(111)
        ax.plot(self.time/1000, self.Is)
        plt.xlabel('Time (sec)')
        plt.ylabel('Is (nA)')  
        plt.title('Neuron Signal')
        ax.grid('on')
        # plt.show() # CHJ
        plt.show(block=True)
        

## Association class for synaptic conductance inputs for motoneurons        
class SynConSignalGenerator:
    def __init__(self, MP, unitCount, department, synType='Excitatory', signalType='Ramp', heav_param=[], iValue=0., pValue=0., period=10000, tau=0., std_max=0., noise=False, t_stop=20000, t_dt=0.1):
        self.genSignal(synType, signalType, heav_param, iValue, pValue, period, tau, std_max, noise, t_stop, t_dt)
        iValue, pValue = pValue, iValue
        synType='Inhibitory'
        self.genSignal(synType, signalType, heav_param, iValue, pValue, period, tau, std_max, noise, t_stop, t_dt)
        self.setSignal(MP, 1, unitCount, department)
        
    def heav(self, x):
        return np.where(x < 0, 0, 1)
        
    def importSignalFile(self, synType, filePath):
        data=pd.read_csv(filePath)
        
        if(synType=='Excitatory'):
            self.G_e = np.asarray(data.G)
            self.e_times = np.asarray(data.time)
        elif(synType=='Inhibitory'):
            self.G_i = np.asarray(data.G)
            self.i_times = np.asarray(data.time)
        
    def genSignal(self, synType, signalType, heav_param, iValue, pValue, period, tau, std_max, noise, t_stop, dt):
        num_steps = int(np.floor(t_stop/dt+1))
        t = np.linspace(0, t_stop, num_steps)

        # excitatory synaptic input
        if(synType=='Excitatory'):
            G_e = np.zeros(num_steps)            
            self.e_times = t    

            if(noise == False):
                if(signalType=='Step'):    
                    i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = heav_param

                    G_e = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                            + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                            + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                            + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                            + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
                    
                elif(signalType=='Ramp'):
                    iv=float(iValue)
                    pv=float(pValue)
                    p=float(period)
                    G_e = pv-((pv-iv)/p)*np.abs(t-p)

            else:
                G_e1 = np.zeros(num_steps)
                amp_e = np.zeros(num_steps)
                tau_e = tau # 0.5
                std_e = std_max # 0.03
                
                if(tau_e !=0):
                    if(signalType=='Step'):
                        i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = heav_param
                        G_e0 = max([i0, ip1, ip2, ip3, ip4, ip5])
                        norm_e = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                                + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                                + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                                + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                                + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
                        norm_e /= G_e0 # normalization
                        
                    elif(signalType=='Ramp'):
                        norm_e = np.zeros(num_steps) # normalization
                        iv = float(iValue)
                        pv = float(pValue)
                        p = float(period)
                        
                        if(pv >= iv):
                            G_e0 = pv # G_max
                            if(pv == 0):
                                norm_e = np.ones(num_steps)
                            else:
                                x = iv/pv
                                norm_e =1-((1-x)/p)*np.abs(t-p) 
                        else:
                            G_e0 = iv # G_max
                            x = pv/iv
                            norm_e =x-((x-1)/p)*np.abs(t-p)
                    
                    # truncated gaussian distribution
                    # for i in xrange(0, num_steps): # CHJ
                    for i in range(0, num_steps):
                        if(norm_e[i] < 0):
                            norm_e[i] = 0.
                                
                    # vectorization
                    exp_e = exp(-dt/tau_e)
                    amp_e = np.sqrt(norm_e)*std_e * np.sqrt( (1-np.exp(-2*dt/tau_e)) )

                    # exact update rule
                    # for i in xrange(1, num_steps): # CHJ
                    for i in range(1, num_steps):
                        G_e1[i] = exp_e * G_e1[i - 1] + amp_e[i - 1] * np.random.normal(loc=0.0, scale=1.0)  # g(t+dt) = g(t) * exp(-dt/tau) + A * N(0,1)
                    G_e = norm_e * G_e0 + G_e1
                
                # white noise
                else:
                    # for i in xrange(0, num_steps): # CHJ
                    for i in range(0, num_steps):
                        G_e[i] = std_e * np.random.normal(loc=0.0, scale=1.0)
            
            # truncated gaussian distribution
            # for i in xrange(0, num_steps): # CHJ
            for i in range(0, num_steps):
                if(G_e[i] < 0):
                    G_e[i] = 0.
                    
            self.G_e = G_e
            
        # inhibitory synaptic input
        if(synType=='Inhibitory'):
            self.i_times = t

            if(noise == False):
                if(signalType=='Step'):
                    i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = heav_param

                    G_i = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                            + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                            + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                            + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                            + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
                    
                elif(signalType=='Ramp'):
                    iv = float(iValue)
                    pv = float(pValue)
                    p = float(period)
                    G_i = pv-((pv-iv)/p)*np.abs(t-p)

            else:
                G_i1 = np.zeros(num_steps)
                amp_i = np.zeros(num_steps)
                tau_i = tau# 2
                std_i = std_max # 0.06
                
                if(tau_i !=0):
                    if(signalType=='Step'):
                        i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = heav_param
                        G_i0 = max([i0, ip1, ip2, ip3, ip4, ip5])
                        norm_i = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                                + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                                + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                                + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                                + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
                        norm_i /= G_i0 # normalization
                        
                    elif(signalType=='Ramp'):
                        norm_i = np.zeros(num_steps) # normalization
                        iv = float(iValue)
                        pv = float(pValue)
                        p = float(period)
                        
                        # vectorization
                        # ramp
                        if(pv >= iv):
                            G_i0 = pv # G_max
                            if(pv == 0):
                                norm_i = np.ones(num_steps)
                            else:
                                x = iv/pv
                                norm_i =1-((1-x)/p)*np.abs(t-p)
                        
                        # inverse ramp
                        else:
                            G_i0 = iv # G_max
                            x = pv/iv
                            norm_i =x-((x-1)/p)*np.abs(t-p) # 1-->x-->1
                            
                    # truncated gaussian distribution
                    # for i in xrange(0, num_steps): # CHJ
                    for i in range(0, num_steps):
                        if(norm_i[i] < 0):
                            norm_i[i] = 0.
                            
                    # vectorization
                    exp_i = exp(-dt/tau_i)
                    amp_i = np.sqrt(norm_i)*std_i * np.sqrt( (1-np.exp(-2*dt/tau_i)) )
                    
                    # exact update rule
                    # for i in xrange(1, num_steps): # CHJ
                    for i in range(1, num_steps):
                        G_i1[i] = exp_i * G_i1[i - 1] + amp_i[i - 1] * np.random.normal(loc=0.0, scale=1.0)
                    G_i = norm_i * G_i0 + G_i1
                
                # white noise
                else:
                    # for i in xrange(0, num_steps): # CHJ
                    for i in range(0, num_steps):
                        G_i[i] = std_i * np.random.normal(loc=0.0, scale=1.0)
            
            # for i in xrange(0, num_steps): # CHJ
            for i in range(0, num_steps):
                if(G_i[i] < 0):
                    G_i[i] = 0.
            
            self.G_i = G_i

    def setSignal(self, MP, MN_range_1, MN_range_2, department):
        MP.setSynConSignals(self.G_e, self.G_i, self.e_times, self.i_times, MN_range_1, MN_range_2, department)

    def plotSignal(self, department):

        fig = plt.figure(dpi=100) 
        # fig.canvas.set_window_title('Synaptic Conductance Signal - ' + department) # CHJ
        # fig.canvas.manager.window.setWindowTitle('Synaptic Conductance Signal - ' + department) # CHJ # for Qt
        fig.canvas.manager.set_window_title('Synaptic Conductance Signal - ' + department) # CHJ # for Inline
        
        ax1=fig.add_subplot(211)
        ax2=fig.add_subplot(212)
        
        # ax1.plot(self.e_times/1000, self.G_e)
        ax1.plot(self.e_times/1000, self.G_e, color='red') #v4.0
        ax1.set_title('Excitatory')
        
        # ax2.plot(self.i_times/1000, self.G_i)
        ax2.plot(self.i_times/1000, self.G_i, color='blue') #v4.0
        ax2.set_title('Inhibitory')
        
        xlabel = 'Time (sec)'
        ylabel = 'G (mS/cm^2)'

        ax2.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)
        ax2.set_ylabel(ylabel)
        ax1.grid('on')
        ax2.grid('on')

        # plt.show() # CHJ
        plt.show(block=True)


## Association class for length variation of muscle-tendon unit        
class MuscleLengthSignalGenerator:
    def __init__(self, MP, unitCount, signalType='Isometric', ivalue=-8, t_stop=2000., t_dt=0.1):
        self.signalType = signalType
        self.genSignal(signalType, 0., t_stop, ivalue, ivalue, t_stop, t_dt)
        self.setSignal(MP, 1, unitCount)

    def importSignalFile(self, filePath):
        self.signalType = 'Exp'
        data = pd.read_csv(filePath)
        self.times = data.time
        self.xm = signal.savgol_filter(np.array(data.Xm), 49, 3)

    def genSignal(self, signalType, itime, ftime, ivalue, fvalue, t_stop, t_dt):
        self.signalType = signalType
        num=int((t_stop)/t_dt)+1
        times=np.linspace(0, t_stop, num)
        self.times=times
        self.xm = np.zeros(num)
        k=0
        
        if(signalType=='Isokinetic'):
            vm=float(fvalue-ivalue)/float(ftime-itime)
            for t in times:
                if(t<=itime):
                    self.xm[k]=ivalue
                elif(t>itime and t<ftime):
                    self.xm[k]=vm*(t-itime)+ivalue
                else:
                    self.xm[k]=fvalue
                k+=1
                
        elif(signalType=='Isometric'):
                self.xm += ivalue
                
        elif(signalType=='Random'):
            # sampling rate
            ms=0.001
            sample_rate=10000
            t_dt=1./sample_rate

            # filter order and sampling correction
            filter_order = 3000*(t_stop)*ms
            # numtaps = filter_order+1 # CHJ
            numtaps = int(filter_order+1) 
            warmup = filter_order/2
            delay = (warmup) / (sample_rate*ms)
            
            # arrary for sampling times
            t_start=0 # msec
            t_stop=t_stop+(delay) # msec
            num=int(((t_stop - t_start)*ms)/t_dt)+1
            times=np.linspace(t_start,t_stop,num)
            
            # random number generation
            n=25
            low=8*n # mms
            high=-8*n # mm
            y=np.random.uniform(low, high, num)
            
            # FIR filtering
            nyq_rate = sample_rate / 2.
            cutoff_hz = 5
            fir_coeff = signal.firwin(numtaps, cutoff_hz/nyq_rate, window='blackman')
            filtered_signal = signal.lfilter(fir_coeff, 1.0, y)
            xm=filtered_signal[int(warmup):]-8
            times_xm=times[int(warmup):]-delay
            
            self.times = times_xm
            self.xm = xm
            
    def setSignal(self, MP, MF_range_1, MF_range_2):
        MP.setMuscleLengthSignals(self.signalType, self.times, self.xm, MF_range_1, MF_range_2)

    def plotSignal(self):
        fig = plt.figure(dpi=100) 
        ax=fig.add_subplot(111)
        ax.plot(self.times/1000, self.xm)
        plt.xlabel('Time (sec)')
        plt.ylabel('Length (mm)')  
        plt.title('MTU length signal')
        ax.grid('on')
        # plt.show() # CHJ
        plt.show(block=True)


## Association class for current impulse stimulation for muscle-tendon unit        
class SpikeSignalGenerator:  
    def __init__(self, MP, unitCount, signalType='Random', freq=100., scale=0., t_start=100., t_stop=1400., t_dt=0.1):
        self.genSignal(signalType, t_start, t_stop, freq, scale, t_dt)
        self.setSignal(MP, 1, unitCount)

    def importSignalFile(self, filePath, t_dt):
        df_sp=pd.read_csv(filePath)
        self.SpikeTimes = np.array(df_sp.time)
        self.t_dt = t_dt
        self.t_stop = self.SpikeTimes[-1]
        
    def genSignal(self, signalType, t_start, t_stop, freq, scale, t_dt):
        ms=1000
        self.t_stop = t_stop
        self.t_dt = t_dt
                
        if(signalType=='Random'):
            num_steps = int(np.floor((t_stop-t_start)/t_dt)+1)  
            t = np.linspace(0, t_stop, num_steps)
            
            if(freq == 0):
                self.SpikeTimes = np.array([])
                return
                
            period = 1./freq*ms
            
            times = np.arange(t_start, t_stop, period)
            
            if(t_stop == times[-1]+period):
                times = np.append(times, t_stop)
                
            # generation of random stimulation times
            num = len(times)
            rn = np.zeros(num)
            k = 0
            
            if(scale != 0):
                for t in times:
                    if(k == 0):
                        rn[k] = times[k]
                        
                    else:
                        rn[k] = np.random.normal(t, scale, 1)
                    k += 1
                        
            else:
                rn = times
        
            rn = np.around(rn, decimals = 1)
            rn.sort()
            self.SpikeTimes = rn
            
            # fitting to simulation time
            bool_arr = self.SpikeTimes <= self.t_stop
            self.SpikeTimes = self.SpikeTimes[bool_arr]

    def setSignal(self, MP, MF_range_1, MF_range_2):
        MP.setSpikeSignals(self.SpikeTimes, MF_range_1, MF_range_2)

    def plotSignal(self):
        num=int((self.t_stop/self.t_dt)+1)
        times=np.linspace(0, self.t_stop, num)

        x_data=times
        y_data=np.zeros(num)

        for i in self.SpikeTimes:
            k=int((i/self.t_dt))
            y_data[k]=1
                
        fig = plt.figure(dpi=100) 
        ax=fig.add_subplot(111)
        ax.plot(x_data/1000, y_data)
        plt.xlabel('Time (sec)')
        plt.ylabel('Current (a.u.)')  
        plt.title('Impulse current Signal')
        ax.grid('on')
 
        # plt.show() # CHJ
        plt.show(block=True)
        

## Independent class for generation of model parameter values for heterogeneous populations
class ParametersGenerator:
    def __init__(self):
        # range parameters for motoneuron
        self.RN = []
        self.Dpath = []
        self.const_Dpath = False
        self.gms = []
        self.gmd = []
        self.gc = []
        self.cmd = []
        self.cms = []
        self.sf = []
        self.sgna = []
        self.dgcal = []
        self.SNM = []
        self.dgkca = []
        
        # range parameters for muscle-tendon unit
        self.p0 = [] 
        self.tau1 = [] 
        self.tau2 = [] 
        self.KSE = []
        self.AM = []
        self.LM = []
        self.cv = []
        self.phi1 = []
        self.phi3 = []
        self.C1i = []
        self.C1n1 = []
        self.C1n4 = []
        self.C2i = []
        self.C2n1 = []
        self.C2n4 = []
        self.C3 = []
        self.C4 = []
        self.C5 = []
        self.alpha_i = []
        self.beta = []
        self.gamma = []
        self.g1 = []
        self.g2 = []
        self.a0 = []        
        self.b0 = []
        self.c0 = []        
        self.d0 = []
    
    def importParameters(self, what, filepath):
        y_arr = pd.read_csv(filepath)
        
        if(what == 'RN'):
            y_arr = np.asarray(y_arr.RN)
            self.RN = y_arr
        elif(what == 'Dpath'):
            y_arr = np.asarray(y_arr.Dpath)
            self.Dpath = y_arr
        # 241211
        elif(what == 'gms'):
            y_arr = np.asarray(y_arr.gms)
            self.gms = y_arr
        elif(what == 'gmd'):
            y_arr = np.asarray(y_arr.gmd)
            self.gmd = y_arr
        elif(what == 'gc'):
            y_arr = np.asarray(y_arr.gc)
            self.gc = y_arr
        elif(what == 'cmd'):
            y_arr = np.asarray(y_arr.cmd)
            self.cmd = y_arr
        elif(what == 'cms'):
            y_arr = np.asarray(y_arr.cms)
            self.cms = y_arr
        elif(what == 'sf'):
            y_arr = np.asarray(y_arr.sf)
            self.sf = y_arr
        elif(what == 'sgna'):
            y_arr = np.asarray(y_arr.sgna)
            self.sgna = y_arr
        elif(what == 'dgcal'):
            y_arr = np.asarray(y_arr.dgcal)
            self.dgcal = y_arr
        elif(what == 'dgkca'):
            y_arr = np.asarray(y_arr.dgkca)
            self.dgkca = y_arr
        elif(what == 'SNM'):
            y_arr = np.asarray(y_arr.SNM)
            self.SNM = y_arr
        # 241211
        elif(what == 'p0'):
            y_arr = np.asarray(y_arr.p0)
            self.p0 = y_arr
        elif(what == 'tau1'):
            y_arr = np.asarray(y_arr.tau1)
            self.tau1 = y_arr
        elif (what == 'tau2'):
            y_arr = np.asarray(y_arr.tau2)
            self.tau2 = y_arr
        elif (what == 'KSE'):
            y_arr = np.asarray(y_arr.KSE)
            self.KSE = y_arr
        elif (what == 'AM'):
            y_arr = np.asarray(y_arr.AM)
            self.AM = y_arr
        elif (what == 'LM'):
            y_arr = np.asarray(y_arr.LM)
            self.LM = y_arr
        elif (what == 'cv'):
            y_arr = np.asarray(y_arr.cv)
            self.cv = y_arr
        elif (what == 'phi1'):
            y_arr = np.asarray(y_arr.phi1)
            self.phi1 = y_arr
        elif (what == 'phi3'):
            y_arr = np.asarray(y_arr.phi3)
            self.phi3 = y_arr
        elif (what == 'C1i'):
            y_arr = np.asarray(y_arr.C1i)
            self.C1i = y_arr
        elif (what == 'C1n1'):
            y_arr = np.asarray(y_arr.C1n1)
            self.C1n1 = y_arr
        elif (what == 'C1n4'):
            y_arr = np.asarray(y_arr.C1n4)
            self.C1n4 = y_arr
        elif (what == 'C2i'):
            y_arr = np.asarray(y_arr.C2i)
            self.C2i = y_arr
        elif (what == 'C2n1'):
            y_arr = np.asarray(y_arr.C2n1)
            self.C2n1 = y_arr
        elif (what == 'C2n4'):
            y_arr = np.asarray(y_arr.C2n4)
            self.C2n4 = y_arr
        elif (what == 'C3'):
            y_arr = np.asarray(y_arr.C3)
            self.C3 = y_arr
        elif (what == 'C4'):
            y_arr = np.asarray(y_arr.C4)
            self.C4 = y_arr
        elif (what == 'C5'):
            y_arr = np.asarray(y_arr.C5)
            self.C5 = y_arr
        elif (what == 'alpha_i'):
            y_arr = np.asarray(y_arr.alpha_i)
            self.alpha_i = y_arr
        elif (what == 'beta'):
            y_arr = np.asarray(y_arr.beta)
            self.beta = y_arr
        elif (what == 'gamma'):
            y_arr = np.asarray(y_arr.gamma)
            self.gamma = y_arr
        elif (what == 'g1'):
            y_arr = np.asarray(y_arr.g1)
            self.g1 = y_arr
        elif (what == 'g2'):
            y_arr = np.asarray(y_arr.g2)
            self.g2 = y_arr
        elif (what == 'a0'):
            y_arr = np.asarray(y_arr.a0)
            self.a0 = y_arr
        elif (what == 'b0'):
            y_arr = np.asarray(y_arr.b0)
            self.b0 = y_arr
        elif (what == 'c0'):
            y_arr = np.asarray(y_arr.c0)
            self.c0 = y_arr
        elif (what == 'd0'):
            y_arr = np.asarray(y_arr.d0)
            self.d0 = y_arr
        else:
            y_arr = np.asarray(y_arr.Value)
        
        return y_arr
    
    # 241210
    #def genParameters(self, what, func, x_arr, params, MNs=0):
    def genParameters(self, what, func, x_arr, params, MNs=0, y_range=None):
        if(y_range == None):
            
            # exponential functions
            if(what == 'RN'):
                y_range = [0.4, 4.0] # (fast, slow)
            elif(what == 'Dpath'):
                y_range =[0.1, 1.1] # (fast, slow)
            # 241211
            elif(what == 'gms'):
                y_range =[0.005, 0.762] # (fast, slow)
            elif(what == 'gmd'):
                y_range =[0.039, 1.840] # (slow, fast)
            elif(what == 'gc'):
                y_range =[0.010, 1.242] # (fast, slow)
            elif(what == 'cms'):
                y_range =[0.125, 4.422] # (slow, fast)
            elif(what == 'cmd'):
                y_range =[0.234, 1.136] # (fast, slow)
            elif(what == 'sf'):
                y_range =[0.005, 0.2] # (fast, slow)
            elif(what == 'sgna'):
                y_range =[30.7, 144.5] # (slow, fast)
            elif(what == 'dgcal'):
                y_range =[0.104, 3.150] # (fast, slow)
            elif(what == 'dgkca'):
                y_range =[0.001, 0.150] # (slow, fast)
            elif(what == 'SNM'):
                y_range =[0.1, 10] # (fast, slow)
            # 241211
            elif(what == 'p0'):
                y_range = [0.1, 10.0] # (slow, fast)
            elif(what == 'tau1'):
                y_range = [1.0, 3.0] # (fast, slow)
            elif(what == 'tau2'):
                y_range = [13.0, 25.0] # (fast, slow)
            elif(what == 'KSE'):
                y_range = [0.16, 0.4] # (fast, slow)
            elif(what == 'AM'):
                y_range = [0.1, 0.5] # (slow, fast)
            elif(what == 'LM'):
                y_range = [1.0, 0.5] # (slow, fast)
            elif(what == 'cv'):
                # 260211                
                y_range = [57.0, 117.0] # (slow, fast)
                # 260211
            elif (what == 'phi1'):
                y_range = [0.002, 0.03] # (fast, slow)
            elif (what == 'phi3'):
                y_range = [1.e-4, 0.01] # (fast, slow)
            elif (what == 'C1i'):
                y_range = [0.12, 0.16] # (fast, slow)
            elif (what == 'C1n1'):
                y_range = [0, 0.01] # (slow, fast)
            elif (what == 'C1n4'):
                y_range = [1.e-4, 85] # (slow, fast)
            elif (what == 'C2i'):
                y_range = [0.09, 0.15]  # (slow, fast)
            elif (what == 'C2n1'):
                y_range = [-0.04, 0.]  # (fast, slow)
            elif (what == 'C2n4'):
                y_range = [1.e-4, 70]  # (slow, fast)
            elif (what == 'C3'):
                y_range = [54, 62]  # (fast, slow)
            elif (what == 'C4'):
                y_range = [-19, -13]  # (fast, slow)
            elif (what == 'C5'):
                y_range = [3.9, 5.1]  # (fast, slow)
            elif (what == 'alpha_i'):
                y_range = [1.6, 2]  # (fast, slow)
            elif (what == 'beta'):
                y_range = [0.09, 0.5]  # (fast, slow)
            elif (what == 'gamma'):
                y_range = [0.0002, 0.001]  # (fast, slow)
            elif (what == 'g1'):
                y_range = [-8, -0.8]  # (slow, fast)
            elif (what == 'g2'):
                y_range = [17, 22]  # (fast, slow)
            elif (what == 'a0'):
                y_range = [0.004, 0.1]  # (fast, slow)
            elif (what == 'b0'):
                y_range = [24, 100]  # (slow, fast)
            elif (what == 'c0'):
                y_range = [-0.58, -0.32]  # (fast, slow)
            elif (what == 'd0'):
                y_range = [30, 43]  # (slow, fast)

        # linear function
        if(func == 'linear'):
            slope, offset = params
    
            y_arr = slope * x_arr + offset
            set_y = set(y_arr)
            # print y_arr # CHJ
            print(y_arr)
            if(len(set_y) == 1):
                self.const_Dpath = True
            else:
                self.const_Dpath = False

        # nonlinear (exponential) functions
        else:
            domain_1, domain_2 = params
            exp_base = np.e
            # 241210            
            #x_arr = np.linspace(domain_1, domain_2, num=MNs)
            x_temp = np.linspace(domain_1, domain_2, num=MNs)

            # increasing convex function
            if(func == 'inc_convex'):
                # 241210                
                #y_arr = exp_base ** x_arr
                y_arr = exp_base ** x_temp
                zero_offset = y_arr[0]
                y_arr = y_arr - zero_offset
                zero_offset = y_range[0]
                gain = abs(y_range[1] - zero_offset) / y_arr[-1]
                y_arr = gain * y_arr + zero_offset
                # 241210                
                spline = splrep(np.linspace(min(x_arr), max(x_arr), num=MNs), y_arr, k = 1)
                y_arr = splev(x_arr, spline)

            # increasing concave function
            if(func == 'inc_concave'):
                # 241210                
                #y_arr = -exp_base ** -x_arr
                y_arr = -exp_base ** -x_temp
                zero_offset = y_arr[0]
                y_arr = y_arr - zero_offset
                zero_offset = y_range[0]
                gain = abs(y_range[1] - zero_offset) / y_arr[-1]
                y_arr = gain * y_arr + zero_offset
                # 241210                
                spline = splrep(np.linspace(min(x_arr), max(x_arr), num=MNs), y_arr, k = 1)
                y_arr = splev(x_arr, spline)

            # decreasing convex function
            if(func == 'dec_convex'):
                # 241210                
                #y_arr = exp_base ** -x_arr
                y_arr = exp_base ** -x_temp
                zero_offset = -y_arr[-1]
                y_arr = y_arr + zero_offset
                zero_offset = y_range[0]
                gain = abs(y_range[1] - zero_offset) / y_arr[0]
                y_arr = gain * y_arr + zero_offset
                # 241210                
                spline = splrep(np.linspace(min(x_arr), max(x_arr), num=MNs), y_arr, k = 1)
                y_arr = splev(x_arr, spline)

            # decreasing concave function
            if(func == 'dec_concave'):
                # 241210                
                #y_arr = -exp_base ** x_arr
                y_arr = -exp_base ** x_temp
                zero_offset = -y_arr[-1]
                y_arr = y_arr + zero_offset
                zero_offset = y_range[0]
                gain = abs(y_range[1]-zero_offset) / y_arr[0]
                y_arr = gain * y_arr + zero_offset
                # 241210                
                spline = splrep(np.linspace(min(x_arr), max(x_arr), num=MNs), y_arr, k = 1)
                y_arr = splev(x_arr, spline)
            
        if(what == 'RN'):
            self.RN = y_arr
        elif(what == 'Dpath'):
            self.Dpath = y_arr
        # 241211
        elif(what == 'gms'):
            self.gms = y_arr
        elif(what == 'gmd'):
            self.gmd = y_arr
        elif(what == 'gc'):
            self.gc = y_arr
        elif(what == 'cms'):
            self.cms = y_arr
        elif(what == 'cmd'):
            self.cmd = y_arr
        elif(what == 'sf'):
            self.sf = y_arr
        elif(what == 'sgna'):
            self.sgna = y_arr
        elif(what == 'dgcal'):
            self.dgcal = y_arr
        elif(what == 'dgkca'):
            self.dgkca = y_arr
        elif(what == 'SNM'):
            self.SNM = y_arr
        # 241211
        elif(what == 'p0'):
            self.p0 = y_arr
        elif(what == 'tau1'):
            self.tau1 = y_arr
        elif (what == 'tau2'):
            self.tau2 = y_arr
        elif (what == 'KSE'):
            self.KSE = y_arr
        elif (what == 'AM'):
            self.AM = y_arr
        elif (what == 'LM'):
            self.LM = y_arr
        elif (what == 'cv'):
            self.cv = y_arr
        elif (what == 'phi1'):
            self.phi1 = y_arr
        elif (what == 'phi3'):
            self.phi3 = y_arr
        elif (what == 'C1i'):
            self.C1i = y_arr
        elif (what == 'C1n1'):
            self.C1n1 = y_arr
        elif (what == 'C1n4'):
            self.C1n4 = y_arr
        elif (what == 'C2i'):
            self.C2i = y_arr
        elif (what == 'C2n1'):
            self.C2n1 = y_arr
        elif (what == 'C2n4'):
            self.C2n4 = y_arr
        elif (what == 'C3'):
            self.C3 = y_arr
        elif (what == 'C4'):
            self.C4 = y_arr
        elif (what == 'C5'):
            self.C5 = y_arr
        elif (what == 'alpha_i'):
            self.alpha_i = y_arr
        elif (what == 'beta'):
            self.beta = y_arr
        elif (what == 'gamma'):
            self.gamma = y_arr
        elif (what == 'g1'):
            self.g1 = y_arr
        elif (what == 'g2'):
            self.g2 = y_arr
        elif (what == 'a0'):
            self.a0 = y_arr
        elif (what == 'b0'):
            self.b0 = y_arr
        elif (what == 'c0'):
            self.c0 = y_arr
        elif (what == 'd0'):
            self.d0 = y_arr
            
        return y_arr
        
    def calCableParameters(self, RN, p, Dpath=0, tm=None, VA=[]):
        if(tm == None):
            # membrane time constant 
            tm = 2.9*RN + 4.5
            tm *= 1.e-3 # unit conversion from ms to s
        
            # voltage attenuation factors
            lamda_sd_dc = -9.32*RN**3 + 82.36*RN**2 + 680.18*RN + 6.82
            lamda_sd_dc *= 1.e-3 # unit conversion from um to mm
            
            lamda_sd_ac = 13.49*RN**3 - 110.43*RN**2 + 286.47*RN + 188.03
            lamda_sd_ac *= 1.e-3 # unit conversion from um to mm
            
            alpha1 = 44.74*RN**3 - 369.33*RN**2 + 1050*RN + 118.35
            alpha1 *= 1.e-3 # unit conversion from um to mm
            
            alpha2 = 8.71*RN**3 - 69.68*RN**2 + 167.94*RN + 111.56
            alpha2 *= 1.e-3 # unit conversion from um to mm
            
            va_sd_dc = np.exp(-Dpath/lamda_sd_dc)
            va_sd_ac = np.exp(-Dpath/lamda_sd_ac)
            va_ds_dc = 1/( 1-np.exp(-alpha1/alpha2)+np.exp((Dpath-alpha1)/alpha2) )

            va_sd_dc = np.round(va_sd_dc, 3)
            va_sd_ac = np.round(va_sd_ac, 3)
            va_ds_dc = np.round(va_ds_dc, 3)

        # initialization of voltage attenuation factors
        else:
            va_sd_dc, va_sd_ac, va_ds_dc = VA
                
        rn = RN * 0.31576 # conversion of input resistance (RN, MOhm) to specific input resistance (rn, ohm * m^2) 
        
        w = 1570.796 # frequency of a single action potential (2*pie*(1/0.004), rad/s)
        
        gms=(1-va_ds_dc)/(rn*(1-va_sd_dc*va_ds_dc))
        
        gmd=(p*va_ds_dc*(1.-va_sd_dc))/((1.-p)*rn*va_sd_dc*(1.-va_sd_dc*va_ds_dc))
    
        gc=(p*va_ds_dc)/(rn*(1.-va_sd_dc*va_ds_dc))

        cmd=(1./(w*(1-p)))*np.sqrt(((gc**2)/(va_sd_ac**2))-((gc+gmd*(1-p))**2))
        
        cms=(tm*(p*(1-p)*tm*gms*gmd+p*gms*(tm*gc-cmd)+(p**2)*gms*cmd+(1-p)*(tm*gc*gmd-gc*cmd)))/(p*((1-p)*(tm*gmd-cmd)+(tm*gc)))
        
        gms = gms*(1.e-1) # unit conversion from S/m^2 to mS/cm^2
        gmd = gmd*(1.e-1) # unit conversion from S/m^2 to mS/cm^2
        gc = gc*(1.e-1) # unit conversion from S/m^2 to mS/cm^2
        cmd = cmd*(1.e+2) # unit conversion from F/m^2 to uF/cm^2
        cms = cms*(1.e+2) # unit conversion from F/m^2 to uF/cm^2
        
        self.gms = gms
        self.gmd = gmd
        self.gc = gc
        self.cms = cms
        self.cmd = cmd
        
        return gms,gmd,gc,cmd,cms
    
    def interpParameters(self, unitType, RN, Dpath=None):
        RN_org = [0.4, 0.5, 0.65, 0.85, 1.1, 1.4, 1.9, 2.6, 3.3, 4]
        
        if(unitType=='motoneuron'):
            # lookup table with linear relationship: y = f(RN_org)
            # sf & sgna
            sf_y = [0.2, 0.19, 0.18, 0.175, 0.015, 0.01, 0.009, 0.007, 0.006, 0.005]
            sgna_y = [54.4, 41.5, 38.7, 47.3, 71.6, 112.2, 144.5, 105.9, 58., 32.6]
            table_y_2Darr = [sf_y, sgna_y]
            # linear interpolation
            k=1
            table_x = RN_org
            interp_x = RN
            i = 0
            for table_y in table_y_2Darr:
                table_y_2Darr[i] = self.cal_interp(table_x, table_y, k, interp_x)
                i += 1
            
            # dgcal
            x = np.array(RN_org[::-1]) # reverse RN array
            y = np.array([0.1, 0.3 ,0.5, 0.6, 0.7, 0.9, 1.1]) # Dpath
            # MN           1      2      3      4      5      6      7      8      9      10        Dpath
            # RN           4.0    3.3    2.6    1.9    1.4    1.1    0.85   0.65   0.5    0.4          
            Z = np.array([[0.104, 0.104, 0.105, 0.107, 0.109, 0.112, 0.118, 0.128, 0.143, 0.165],  #  0.1
                          [0.107, 0.110, 0.113, 0.119, 0.128, 0.141, 0.165, 0.220, 0.340, 0.578],  #  0.3
                          [0.110, 0.112, 0.119, 0.132, 0.147, 0.172, 0.218, 0.335, 0.670, 1.340],  #  0.5
                          [0.111, 0.113, 0.120, 0.134, 0.157, 0.185, 0.241, 0.396, 0.852, 1.665],  #  0.6
                          [0.112, 0.117, 0.126, 0.144, 0.172, 0.206, 0.269, 0.449, 0.958, 1.815],  #  0.7
                          [0.123, 0.135, 0.147, 0.168, 0.212, 0.276, 0.416, 0.770, 1.510, 2.110],  #  0.9
                          [0.182, 0.218, 0.238, 0.270, 0.389, 0.535, 0.550, 0.767, 1.850, 3.150]]) #  1.1
            # linear interpolation
            f = interp2d(x, y, Z, kind='linear')
            dgcal = f(RN, Dpath)
            RN_idx = np.argsort(np.argsort(RN))
            Dpath_idx = np.argsort(np.argsort(Dpath))
            dgcal = dgcal[:, RN_idx][Dpath_idx, :]
            
            if(self.const_Dpath == True):
                dgcal = dgcal[0,:]
            else:
                dgcal = np.diag(dgcal)
            
            # SNM
            # uniform distribution of brainstem neuromodulation
            SNM_y = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) 
            # linear interpolation
            SNM = self.cal_interp(table_x, SNM_y, k, interp_x)
            dgcal = SNM * dgcal
            
            # dgkca
            x = np.array(RN_org[::-1]) # reverse RN array
            y = np.array([0.1, 0.3 ,0.5, 0.6, 0.7, 0.9, 1.1])
            # uniform distribution of dendritic K(Ca) channels
            # MN           1        2        3        4        5        6        7        8        9        10          Dpath
            # RN           4.0      3.3      2.6      1.9      1.4      1.1      0.85     0.65     0.5      0.4          
            Z = np.array([[0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001],  #  0.1
                          [0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001],  #  0.3
                          [0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001],  #  0.5
                          [0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001],  #  0.6
                          [0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001],  #  0.7
                          [0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001],  #  0.9
                          [0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001]]) #  1.1

            # MN           1        2        3        4        5        6        7        8        9        10          Dpath
            # RN           4.0      3.3      2.6      1.9      1.4      1.1      0.85     0.65     0.5      0.4          
#            Z = np.array([[0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15],  #  0.1
#                          [0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15],  #  0.3
#                          [0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15],  #  0.5
#                          [0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15],  #  0.6
#                          [0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15],  #  0.7
#                          [0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15],  #  0.9
#                          [0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15,    0.15]]) #  1.1
            # linear interpolation
            f = interp2d(x, y, Z, kind='linear')
            dgkca = f(RN, Dpath)
            RN_idx = np.argsort(np.argsort(RN))
            Dpath_idx = np.argsort(np.argsort(Dpath))
            dgkca = dgkca[:, RN_idx][Dpath_idx, :]
            
            if(self.const_Dpath == True):
                dgkca = dgkca[0,:]
            else:
                dgkca = np.diag(dgkca)
            
            # display of range parameters for motoneuron
            self.sf = table_y_2Darr[0]
            self.sgna = table_y_2Darr[1]
            self.dgcal = dgcal
            self.SNM = SNM
            self.dgkca = dgkca
            
            return [i for i in table_y_2Darr] + [dgcal, dgkca]
        
        if(unitType=='musclefibers' or unitType=='motorunit'):
            p0_org = [0.1000, 2.0250, 3.9500, 5.8750, 7.2500, 8.0750, 8.7625, 9.3125, 9.7250, 10]
            
            # 241208            
            # lookup table with linear relationship: p0_org = f(RN_org) for motor unit
            #x = RN_org
            #p0_y = [10, 9.7250, 9.3125, 8.7625, 8.0750, 7.2500, 5.8750, 3.9500, 2.0250, 0.1000]
            #spline = splrep(x, p0_y, k = 1)
            #p0_interp = splev(RN, spline)
            #self.p0 = p0_interp
            
            p0_interp = RN
            self.p0 = p0_interp

            
            # lookup table with linear relationship: y = f(p0_org)
            # tau1 (1.0 ~ 3.0) (fast to slow)
            tau1_y = [3, 2.6111, 2.2222, 1.8333, 1.5556, 1.3889, 1.2500, 1.1389, 1.0556, 1]
            # tau2 (13 ~ 25) (fast to slow)
            tau2_y = [25, 22.6667, 20.3333, 18, 16.3333, 15.3333, 14.5000, 13.8333, 13.3333, 13]
            # KSE (0.16 ~ 0.4) (fast to slow)
            KSE_y = [0.4000, 0.3533, 0.3067, 0.2600, 0.2267, 0.2067, 0.1900, 0.1767, 0.1667, 0.1600]
            # AM (0.5 - 0.1) (fast to slow)
            AM_y = [0.1000, 0.1778, 0.2556, 0.3333, 0.3889, 0.4222, 0.4500, 0.4722, 0.4889, 0.5000]
            # LM (0.5 ~ 1.0) (fast to slow)
            LM_y = [1, 0.9028, 0.8056, 0.7083, 0.6389, 0.5972, 0.5625, 0.5347, 0.5139, 0.5000]
            # 260211            
            # cv (117 - 57) (fast to slow)
            cv_y = [57, 68.67, 80.33, 92.0, 100.33, 105.33, 109.5, 112.83, 115.33, 117.0]            
            # 260211            
            # phi1 (0.002 ~ 0.03) (fast to slow)
            phi1_y = [0.0300, 0.0246, 0.0191, 0.0137, 0.0098, 0.0074, 0.0055, 0.0039, 0.0028, 0.0020]
            # phi3 (0.0001 ~ 0.01) (fast to slow)
            phi3_y = [0.0100, 0.00808, 0.00615, 0.00423, 0.00285, 0.00203, 0.00134, 0.00079, 0.00038, 0.0001]
            # C1i (0.16 ~ 0.12) (fast to slow)
            C1i_y = [0.1200, 0.1278, 0.1356, 0.1433, 0.1489, 0.1522, 0.1550, 0.1572, 0.1589, 0.1600]
            # C1n1 (0.01 ~ 0) (fast to slow)
            C1n1_y = [0, 0.00194, 0.00389, 0.00583, 0.00722, 0.00806, 0.00875, 0.00931, 0.00972, 0.01000]
            # C1n4 (85 ~ 0.0001) (fast to slow)
            C1n4_y = [0.0001, 16.5279, 33.0556, 49.5834, 61.3889, 68.4722, 74.3750,	79.0972, 82.6389, 85]
            # C2i (0.15 ~ 0.09) (fast to slow)
            C2i_y = [0.0900, 0.1017, 0.1133, 0.1250, 0.1333, 0.1383, 0.1425, 0.1458, 0.1483, 0.1500]
            # C2n1 (-0.04 ~ 0) (fast to slow)
            C2n1_y = [0, -0.0078, -0.0156, -0.0233, -0.0289, -0.0322, -0.0350, -0.0372, -0.0389, -0.0400]
            # C2n4 (70 ~ 0.0001) (fast to slow)
            C2n4_y = [0.0001, 13.6112, 27.2223, 40.8334, 50.5556, 56.389, 61.2500, 65.1389, 68.0556, 70.0000]
            # C3 (54 ~ 62) (fast to slow)
            C3_y = [62, 60.4444, 58.8889, 57.3333, 56.2222, 55.5556, 55, 54.5556, 54.2222, 54]
            # C4 (-19 ~ -13) (fast to slow)
            C4_y = [-13.0000, -14.1667, -15.3333, -16.5000, -17.3333, -17.8333, -18.2500, -18.5833, -18.8333, -19]
            # C5 (3.9 ~ 5.1) (fast to slow)
            C5_y = [5.1000, 4.8667, 4.6333, 4.4000, 4.2333, 4.1333, 4.0500, 3.9833, 3.9333, 3.9000]
            # alpha_i (1.6 ~ 2)
            alpha_i_y = [2, 1.9222, 1.8444, 1.7667, 1.7111, 1.6778, 1.6500, 1.6278, 1.6111, 1.6000]
            # beta (0.09 ~ 0.5) (fast to slow)
            beta_y = [0.5000, 0.4203, 0.3406, 0.2608, 0.2039, 0.1697, 0.1413, 0.1185, 0.1014, 0.0900]
            # gamma (0.0002 ~ 0.001) (fast to slow)
            gamma_y = [0.00100, 0.00084, 0.00069, 0.00053, 0.00042, 0.00036, 0.00030, 0.00026, 0.00022, 0.00020]
            # g1 (-0.8 ~ -8) (fast to slow)
            g1_y = [-8, -6.6000, -5.2000, -3.8000, -2.8000, -2.2000, -1.7000, -1.3000, -1.0000, -0.8000]
            # g2 (17 ~ 22) (fast to slow)
            g2_y = [22, 21.0278, 20.0556, 19.0833, 18.3889, 17.9722, 17.6250, 17.3472, 17.1389, 17]
            # a0 (0.004 ~ 0.1) (fast to slow)
            a0_y = [0.1000, 0.0813, 0.0627, 0.0440, 0.0307, 0.0227, 0.0160, 0.0107, 0.0067, 0.0040]
            # b0 (100 ~ 24) (fast to slow)
            b0_y = [24, 38.7778, 53.5556, 68.3333, 78.8889, 85.2222, 90.5000, 94.7222, 97.889, 100]
            # c0 (-0.58 ~ -0.32) (fast to slow)
            c0_y = [-0.3200, -0.3706, -0.4211, -0.4717, -0.5078, -0.5294, -0.5475, -0.5619, -0.5728, -0.5800]
            # d0 (43 ~ 30) (fast to slow)
            d0_y = [30.0000, 32.5278, 35.0556, 37.5833, 39.3888, 40.4722, 41.3750, 42.0972, 42.6389, 43]
            
            table_y_2Darr = [tau1_y,
                             tau2_y,
                             KSE_y,
                             AM_y,
                             LM_y,
                             cv_y,
                             phi1_y,
                             phi3_y,
                             C1i_y,
                             C1n1_y,
                             C1n4_y,
                             C2i_y,
                             C2n1_y,
                             C2n4_y,
                             C3_y,
                             C4_y,
                             C5_y,
                             alpha_i_y,
                             beta_y,
                             gamma_y,
                             g1_y,
                             g2_y,
                             a0_y,
                             b0_y,
                             c0_y,
                             d0_y]
            
            # linear interpolation of range parameter values for muscle-tendon unit
            table_x = p0_org
            interp_x = p0_interp
            k = 1
            i = 0
            for table_y in table_y_2Darr:
                table_y_2Darr[i] = self.cal_interp(table_x, table_y, k, interp_x)
                i += 1
                
            # display of range parameters for muscle-tendon unit
            self.tau1 = table_y_2Darr[0]
            self.tau2 = table_y_2Darr[1]
            self.KSE = table_y_2Darr[2]
            self.AM = table_y_2Darr[3]
            self.LM = table_y_2Darr[4]
            self.cv = table_y_2Darr[5]
            self.phi1 = table_y_2Darr[6]
            self.phi3 = table_y_2Darr[7]
            self.C1i = table_y_2Darr[8]
            self.C1n1 = table_y_2Darr[9]
            self.C1n4 = table_y_2Darr[10]
            self.C2i = table_y_2Darr[11]
            self.C2n1 = table_y_2Darr[12]
            self.C2n4 = table_y_2Darr[13]
            self.C3 = table_y_2Darr[14]
            self.C4 = table_y_2Darr[15]
            self.C5 = table_y_2Darr[16]
            self.alpha_i = table_y_2Darr[17]
            self.beta = table_y_2Darr[18]
            self.gamma = table_y_2Darr[19]
            self.g1 = table_y_2Darr[20]
            self.g2 = table_y_2Darr[21]
            self.a0 = table_y_2Darr[22]
            self.b0 = table_y_2Darr[23]
            self.c0 = table_y_2Darr[24]
            self.d0 = table_y_2Darr[25]
            
            # 241208
            #return [p0_interp] + [i for i in table_y_2Darr]
            return [i for i in table_y_2Darr]
        
    def cal_interp(self, table_x, table_y, k, interp_x):
        spline = splrep(table_x, table_y, k=1)
        y = splev(interp_x, spline)
        return y
    
    def plotParameters(self, MP_arr, what=None, unitType=None):
        if(what == 'RN'):
            fig = plt.figure(dpi=100) 
            ax=fig.add_subplot(111)
            ax.plot(MP_arr, self.RN, '-o')

            if(unitType == 'motorunit'):
                plt.xlabel('MU#')
                plt.ylabel('RN (M*Ohm)')
                plt.title('MU#_RN')
            else:
                plt.xlabel('MN#')
                plt.ylabel('RN (M*Ohm)')
                plt.title('MN#_RN')
            ax.grid('on')
        
        elif(what == 'Dpath'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.Dpath, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('Dpath (mm)')
            plt.title('RN-Dpath')
            ax2.grid('on')
        
        elif(what == 'gms'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.gms, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('gms (mS/cm^2)')
            plt.title('RN-gms')
            ax2.grid('on')
            
        elif(what == 'gmd'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.gmd, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('gmd (mS/cm^2)')
            plt.title('RN-gmd')
            ax2.grid('on')
        
        elif(what == 'gc'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.gc, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('gc (mS/cm^2)')
            plt.title('RN-gc')
            ax2.grid('on')
            
        elif(what == 'cms'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.cms, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('cms (uF/cm^2)')
            plt.title('RN-cms')
            ax2.grid('on')
        
        elif(what == 'cmd'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.cmd, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('cmd (uF/cm^2)')
            plt.title('RN-cmd')
            ax2.grid('on')
            
        elif(what == 'sf'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.sf, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('sf')
            plt.title('RN-sf')
            ax2.grid('on')
        
        elif(what == 'sgna'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.sgna, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('sgna (mS/cm^2)')
            plt.title('RN-sgna')
            ax2.grid('on')
        
        elif(what == 'dgcal'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.dgcal, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('dgcal (mS/cm^2)')
            plt.title('RN-dgcal')
            ax2.grid('on')
        
        elif(what == 'SNM'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.SNM, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('SNM')
            plt.title('RN-SNM')
            ax2.grid('on')
        
        elif(what == 'dgkca'):
            fig = plt.figure(dpi=100) 
            ax2=fig.add_subplot(111)
            ax2.plot(self.RN, self.dgkca, '-o')
            plt.xlabel('RN (M*Ohm)')
            plt.ylabel('dgkca (mS/cm^2)')
            plt.title('RN-dgkca')
            ax2.grid('on')
            
        elif(what == 'p0'):
            fig = plt.figure(dpi=100) 
            ax=fig.add_subplot(111)
            
            if(unitType == 'motorunit'):
                ax.plot(self.RN, self.p0, '-o')
                plt.xlabel('RN (M*Ohm)') #v4.0
                plt.ylabel('P0 (N)')
                plt.title('RN_P0')
                
                # print 'RN: ' # CHJ
                print('RN: ')
                # print self.RN # CHJ
                print(self.RN)
            else:
                ax.plot(MP_arr, self.p0, '-o')
                plt.xlabel('MTU#')
                plt.ylabel('P0 (N)')
                plt.title('MTU#_P0')
            ax.grid('on')

        elif (what == 'tau1'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.tau1, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('tau1 (ms)')
            plt.title('P0_tau1')
            ax.grid('on')

        elif (what == 'tau2'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.tau2, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('tau2 (ms)')
            plt.title('P0_tau2')
            ax.grid('on')

        elif (what == 'KSE'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.KSE, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('KSE (/mm)')
            plt.title('P0_KSE')
            ax.grid('on')

        elif (what == 'AM'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.AM, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('AM (mV)')
            plt.title('P0_AM')
            ax.grid('on')

        elif (what == 'LM'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.LM, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('LM (ms)')
            plt.title('P0_LM')
            ax.grid('on')

        elif (what == 'cv'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.cv, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('cv (m/s)')
            plt.title('P0_cv')
            ax.grid('on')

        elif (what == 'phi1'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.phi1, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('phi1')
            plt.title('P0_phi1')
            ax.grid('on')

        elif (what == 'phi3'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.phi3, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('phi3')
            plt.title('P0_phi3')
            ax.grid('on')

        elif (what == 'C1i'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.C1i, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('C1i')
            plt.title('P0_C1i')
            ax.grid('on')

        elif (what == 'C1n1'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.C1n1, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('C1n1')
            plt.title('P0_C1n1')
            ax.grid('on')

        elif (what == 'C1n4'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.C1n4, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('C1n4')
            plt.title('P0_C1n4')
            ax.grid('on')

        elif (what == 'C2i'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.C2i, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('C2i')
            plt.title('P0_C2i')
            ax.grid('on')

        elif (what == 'C2n1'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.C2n1, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('C2n1')
            plt.title('P0_C2n1')
            ax.grid('on')

        elif (what == 'C2n4'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.C2n4, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('C2n4')
            plt.title('P0_C2n4')
            ax.grid('on')

        elif (what == 'C3'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.C3, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('C3')
            plt.title('P0_C3')
            ax.grid('on')

        elif (what == 'C4'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.C4, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('C4')
            plt.title('P0_C4')
            ax.grid('on')

        elif (what == 'C5'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.C5, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('C5')
            plt.title('P0_C5')
            ax.grid('on')

        elif (what == 'alpha_i'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.alpha_i, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('alpha_i')
            plt.title('P0_alpha_i')
            ax.grid('on')

        elif (what == 'beta'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.beta, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('beta')
            plt.title('P0_beta')
            ax.grid('on')

        elif (what == 'gamma'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.gamma, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('gamma')
            plt.title('P0_gamma')
            ax.grid('on')

        elif (what == 'g1'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.g1, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('g1')
            plt.title('P0_g1')
            ax.grid('on')

        elif (what == 'g2'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.g2, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('g2')
            plt.title('P0_g2')
            ax.grid('on')
        
        elif (what == 'a0'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.a0, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('a0')
            plt.title('P0_a0')
            ax.grid('on')        
    
        elif (what == 'b0'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.b0, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('b0')
            plt.title('P0_b0')
            ax.grid('on')
        
        elif (what == 'c0'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.c0, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('c0')
            plt.title('P0_c0')
            ax.grid('on')

        elif (what == 'd0'):
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(self.p0, self.d0, '-o')
            plt.xlabel('P0 (N)')
            plt.ylabel('d0')
            plt.title('P0_d0')
            ax.grid('on')

        # plt.show() # CHJ
        plt.show(block=True)


## Association class for parallel simulation and data management
class ParallelManager:
    def __init__(self):
        # print "Parallel manager instance constructed." # CHJ
        print("Parallel manager instance constructed.")
     
    def killPPserver(self, ls_server):
        for node in ls_server:
            os.system('ssh -f '+node+' fuser -k -n tcp 60000')

    def startPPserver(self, ls_server):
        curPath=os.path.dirname( os.path.abspath( __file__ ) )+'/'
        fileName='ppserver.sh'
        
        for node in ls_server:
            os.system('ssh -f '+node+' '+curPath+fileName)

    def selectComputeNode(self, ls_Node):
        ppservers = ('*', )
        self.jobServer = pp.Server(ncpus=0, ppservers=ppservers)
        time.sleep(3)
        # ppservers=self.jobServer.get_active_nodes().keys() # CHJ
        ppservers = list(self.jobServer.get_active_nodes().keys())
        # num_node=self.jobServer.get_active_nodes().values() # CHJ
        num_node = list(self.jobServer.get_active_nodes().values())

        for i in range(len(ppservers)):
            if(ppservers[i]=='local'):
                ppservers.pop(i)
                break
        

        for i in range(len(ppservers)):
            ppservers[i]=ppservers[i].replace(":60000", "")

        # initialization and activation of PP servers  
        self.killPPserver(ppservers)
        time.sleep(3)

        self.startPPserver(ls_Node)
        time.sleep(3)

        ppservers= []

        for i in ls_Node:
            ppservers.append(i+':60000')

        # print ppservers # CHJ
        print(ppservers)

    def selectJobServer(self, num=0, ls_server=None):
        ppservers = ('*', ) # all servers available

        if(ls_server==None): # only management node with cores available
            self.jobServer = pp.Server(ncpus=num)

        elif(ls_server=='ALL'): # management node + all computational nodes
            self.jobServer = pp.Server(ncpus=num, ppservers=ppservers)
  
        else: # management node + assigned computational nodes
            self.jobServer = pp.Server(ncpus=num, ppservers=tuple(ls_server))
    
   	    # display of participating node information (IP and node number)
        time.sleep(3)
        # print self.jobServer.get_active_nodes() # CHJ
        print(self.jobServer.get_active_nodes())

    def setJobCount(self, num):
        self.jobCount=num
        
    def submitJob(self, job, arr_units):
        
        results = []
    
        for unit in arr_units:
            j = self.jobServer.submit(job, (unit, ), modules=('time', 'os', 'platform', 'import pandas as pd', 'from scipy import integrate', 'import numpy as np', 'from math import log, exp, pi, asin, sin, tanh, cosh', 'import matplotlib.pyplot as plt',))
            results.append(j) 

        self.jobServer.wait()
        
        result_data=[]
        
        for i in results:
            result_data.append(i())
        
        self.result_data=result_data
                
        # print self.jobServer.print_stats() # CHJ
        print(self.jobServer.print_stats())

    def runSimulation(self, MP):
        arr_units=MP.getUnitsArray()
        self.submitJob(job,arr_units)

    def saveResultData(self, MP, savePath, fileName):
        if not os.path.isdir(savePath):
            os.makedirs(savePath)
        
        for i in self.result_data:
            uniqueNumber=i[0]
            resultArray=i[2]

            fileName_ex = fileName + '_%03d' %uniqueNumber + '.csv'
            
            df = pd.DataFrame(resultArray)
            
            unitType = MP.getUnitType()

            if(unitType == 'motoneuron'): 
                cols = ['Time','G_esyn_dend', 'Firing_rate', 'Is', 'V_soma', '[Ca]_soma','E_Ca_soma','I_Naf_soma','m_Naf_soma','h_Naf_soma','I_Nap_soma','m_Nap_soma','I_Kdr_soma','n_Kdr_soma','I_Kca_soma','I_Can_soma','m_Can_soma','h_Can_soma','I_H_soma','m_H_soma','V_dend','[Ca]_dend','E_Ca_dend','I_Cal_dend','l_Cal_dend','I_Naf_dend','m_Naf_dend','h_Naf_dend','I_Nap_dend','m_Nap_dend','I_Kdr_dend','n_Kdr_dend','I_Kca_dend','m_Kca_dend','I_Can_dend','m_Can_dend','h_Can_dend','I_H_dend','m_H_dend','I_esyn_soma','G_esyn_soma','I_isyn_soma','G_isyn_soma','I_esyn_dend','I_isyn_dend','G_isyn_dend']

            elif(unitType == 'musclefibers'):
                cols = ['Time','A','Am','A_tilde','C1','C2','CaSP','CaSPB','CaSPT','CaSR','CaSRCS','F','MUAP','R','Spike','Vm','XCE','Xm']

            elif(unitType == 'motorunit'):
                cols = ['Time','G_esyn_dend', 'Firing_rate', 'MN_Spike', 'Is', 'V_soma', '[Ca]_soma','E_Ca_soma','I_Naf_soma','m_Naf_soma','h_Naf_soma','I_Nap_soma','m_Nap_soma','I_Kdr_soma','n_Kdr_soma','I_Kca_soma','I_Can_soma','m_Can_soma','h_Can_soma','I_H_soma','m_H_soma','V_dend','[Ca]_dend','E_Ca_dend','I_Cal_dend','l_Cal_dend','I_Naf_dend','m_Naf_dend','h_Naf_dend','I_Nap_dend','m_Nap_dend','I_Kdr_dend','n_Kdr_dend','I_Kca_dend','m_Kca_dend','I_Can_dend','m_Can_dend','h_Can_dend','I_H_dend','m_H_dend','I_esyn_soma','G_esyn_soma','I_isyn_soma','G_isyn_soma','I_esyn_dend','I_isyn_dend','G_isyn_dend','MF_Spike','A','Am','A_tilde','C1','C2','CaSP','CaSPB','CaSPT','CaSR','CaSRCS','F','MUAP','R','Vm','XCE','Xm']

            df_final = df.reindex(columns = cols)
            df_final.to_csv(savePath + fileName_ex) 
            
            # print str(fileName_ex)+" file saved." # CHJ
            print(str(fileName_ex)+" file saved.")

    def plotResultData(self, MP, scope_list):
        result_length=len(self.result_data)
                
        nrows=MP.getUnitCount()
        ncols=1
        gs0 = gridspec.GridSpec(nrows, ncols)

        for scope in scope_list: # scope_list = vs, vd
            k=0

                
            for i in self.result_data:
                cell_result=i[2]
                uniqueNumber=i[0]

                x_data = cell_result['Time']/1000
                y_data=cell_result[scope]
                ax = plt.subplot(gs0[k])

                if(scope=='Firing_rate'):
                    x_data=x_data[cell_result['Firing_rate'] > 0]
                    y_data=cell_result['Firing_rate'][cell_result['Firing_rate'] > 0]
                    plt.plot(x_data, y_data, '.')

                else:
                    plt.plot(x_data, y_data)

                unitType=MP.getUnitType()

                if(unitType == 'motoneuron'):
                    plt.ylabel('MN'+str(uniqueNumber))

                elif(unitType == 'musclefibers'):
                    plt.ylabel('MTU'+str(uniqueNumber))

                elif(unitType == 'motorunit'):
                    plt.ylabel('MU'+str(uniqueNumber))

                ax.set_xlim(0, cell_result['Time'][-1]/1000)
                ax.set_autoscaley_on(True)
                ax.autoscale_view()
                plt.grid('on')
                k+=1
            
            plt.xlabel('Time (sec)')
            figure = plt.gcf()
            # figure.canvas.set_window_title(scope) # CHJ
            figure.canvas.manager.set_window_title(scope)
            # plt.show() # CHJ
            plt.show(block=True)


## Job content
def job(cell):
    cell.getHostname()
    MN_R=cell.solModel()

    return cell.uniqueNumber, cell.hostname, MN_R


## Spike detection for motoneuron
tmp_peak = False
def detect_Spike(SpikeTimes, t, Vs):
    global tmp_peak
    HIGH_THRH_SPIKE = -10 # threshold for spike detection
    S_Detect = False

    if (Vs > HIGH_THRH_SPIKE):
        tmp_peak = True

    elif tmp_peak:
        SpikeTimes.append(t)
        tmp_peak = False
        S_Detect = True

    return SpikeTimes, S_Detect
    
    
## Instantaneous firing rate for motoneuron    
def cal_FiringRate(SpikeTimes):    
    i = len(SpikeTimes)-1
    cal_FiringRate = 0
    
    if(i >= 1): 
        s = 1000
        T1 = SpikeTimes[i-1]
        T = SpikeTimes[i]
        cal_FiringRate = 1/(T-T1)*s
    
    return cal_FiringRate