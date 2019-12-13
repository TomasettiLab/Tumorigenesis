# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 21:45:16 2016

@author: kamel
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 23:22:01 2016

@author: kamel
"""
# Two Carrying capacities one for symmetric division and one for asymetric division
# -*- coding: utf-8 -*-
#cancerMutations2.py
#Feedback Model
#Last Updated: 4/26/2015, Added feature to print out population of individual type 3 clones


import os
from os import path
import argparse
import numpy as np
import numpy.matlib as npmat
import time
import warnings
warnings.filterwarnings("ignore")

class tissue:
    #parameters:
    #N: wild type population size
    #ngrp: number of mutation groups
    #pmute: probability of mutation
    #rho0: slope parameter for asymetric differentiation
    #delta: update time step for simulation
    #pgrp: conditional probability of groups given mutation. If
    #      length=1,all groups are equi-likely 
    #T: maximum lifetime (in weeks)
    #t3Threshold: number of mutations in last type for cancer to
    #             occur.
    #epsilon0: radius of lowest to highest rate intervals for normal
    #          population homeostasis
    #cang: angio-genesis effect (increased resource)
    #capp: apoptosis effect (death rate decrease)
    #cpro: proliferation effect (division increase)
    
    def __init__(self, N = 1e1, dcell = 1e9, ngrp=18, pmute =1e-6,pmuteLOH=1e-4,pmutesuperLOH=10**(-3), rho0 = 1-2e-6, pgrp = [1],
                    tau0=1.0, epsilon0 =1-2e-6, cang = 1., alpha=1., capp=1.5 ,cappInter=1., targ_size=10.,targ_size1h=100., delta=0.5, T=750, t3Threshold = 1,deltaMute=5*1e-6,deltaOnc1=0.001,kap0=1,psym=0.053,pasymdiff=0.88,psymdiff=0.049,deltapsym=10e-3):
        self.Tsind=np.append([0,2,3],np.arange(10,18))
        self.Oncind=np.append([1],np.arange(4,10)) 
        self.PathwaysCF=[{0,1,2},{3}]
        self.PathwaysCS=[{4,5,6},{7},{8,9},{10}]
        self.PathwaysGM=[{11},{12},{13},{14},{15},{16},{17}]                
        self.N = long(N)
        self.crypt_size = 10
        self.nb_crypt = self.N /self.crypt_size 
        self.dcell = long(dcell)
        self.ngrp = ngrp
        #self.pchange = pchange
        self.nmut = 0
        self.rho0 = rho0
        self.tau0 = tau0
        self.alpha = alpha
        self.capp = capp
        self.cappInter=cappInter
        self.sigma_app = float(targ_size)*((self.alpha+1)*targ_size + (self.alpha-1)*self.crypt_size)/((self.alpha-1)*targ_size + (self.alpha+1)*self.crypt_size)
        r = targ_size1h
        self.sigma_app1h = (r**2)/self.crypt_size

        self.cang = cang
        self.epsilon0 = epsilon0
        self.delta = delta
        self.deltapsym=deltapsym
        if len(pgrp) == 1:
            self.pgrp = np.append(np.ones(ngrp-1)/(ngrp-1),0)
            self.pgrp[4] *=1
            self.pgrp =self.pgrp/self.pgrp.sum()
        else:
            self.pgrp = pgrp
            
        self.pmute = pmute
        self.pmuteLOH=pmuteLOH
        self.pmutesuperLOH=pmutesuperLOH
        self.deltaMute=deltaMute
        self.T = T
        self.t3Threshold = t3Threshold
        self.groupsize = np.array([N] + [0]*self.ngrp, dtype=long)
        self.kchange = .1
        self.maxAge = 500
        self.deltaOnc1=deltaOnc1
        self.kap0=kap0
        self.psym=psym
        self.pasymdiff=pasymdiff
        self.psymdiff=psymdiff
        self.Gen_mutations=self.Generate_Mutations()
        self.mutations_prob=self.Generate_Mutations_Prob()
        
        
   

    #basic division/death rate with total population feedback
    #basic division/death rate with total population feedback
    def hfun(self, n, ka, s, ep):
        #a=np.log(1+((s-1)*float(n)/float(s)))/np.log(s)
        #u = float(7)/(3.5+a)
        if (float(10*n)/(2*10**5*s))<1:
           u =2.5*np.exp(1)*0.45/3.5*np.exp(-1.0/(1-float(100*n**2)/(4*10**10*s**2)))+0.45/3.5
        else:
           u=0.45/3.5
        
        
            #a * self.crypt_size + b* float(n))/(self.crypt_size + float(n))
            # s is the target size
            # n is ntot
        return u
    
    
        
    

    #Determine rates and probabilities for a given genotype and
    #possibly age 
    def getProb(self, mutations, age, nc):
        #ntot = float(self.N + self.nmut)
        psym=self.psym
        psymdiff=self.psymdiff
        pasymdiff=self.pasymdiff
        ntot = float(self.crypt_size+nc)
        tauDiv=self.hfun(ntot,1,self.crypt_size,self.epsilon0)
        tauDeathc = self.hfun(self.crypt_size,1,self.crypt_size,self.epsilon0)*(psym-psymdiff)
        pmute=self.pmute
        
        if np.any(mutations[[0,2,3]]>=2) or np.any(mutations[1]>0):
            psym +=self.deltapsym
            psymdiff -=self.deltapsym
            
        
        if np.any(mutations[10]>=2) or np.any(mutations[np.arange(4,10)]>0):
            tauDiv *=self.capp
        
        if np.any(mutations[np.arange(11,18)]>=2):
            pmute=self.deltaMute
        
        
        
        if np.any(mutations[self.Tsind]==1):
            pmute=self.pmuteLOH
        
        tauSum=tauDeathc+tauDiv
        tauNot=max(0,1-tauSum)
        tauDiv /=(tauSum+tauNot)
        tauDeathc /=(tauSum+tauNot)
        tauDeath=tauDeathc+tauDiv*psymdiff
        p=(tauDeath,tauDiv*pasymdiff*(1-float(pmute)/2),tauDiv*pasymdiff*float(pmute)/2,tauDiv*psym*(1-pmute),tauDiv*psym*pmute,tauNot)
        return p
    
    
    def getPreRate(self, mutations, age, nc):
        #ntot = float(self.N + self.nmut)
        mut=mutations+np.append(1,np.zeros(17))
        countCF=0
        countCS=0
        countGM=0
        setTS=(set(np.where(mut>=2)[0]) & set(self.Tsind))
        setOnc=(set(np.where(mut>0)[0]) & set(self.Oncind))
        setDr=setTS or setOnc
        psym=self.psym
        psymdiff=self.psymdiff
        pasymdiff=self.pasymdiff
        ntot = float(self.crypt_size+nc)
        tauDiv=self.hfun(ntot,1,self.crypt_size,self.epsilon0)
        tauDeathc = self.hfun(self.crypt_size,1,self.crypt_size,self.epsilon0)*(psym-psymdiff)
        pmute=(countGM)*self.deltaMute+self.pmute
        
        for i in self.PathwaysCF:
            countCF +=bool(i & setDr)
        
        for i in self.PathwaysCS:
            countCS +=bool(i & setDr)
        
        for i in self.PathwaysGM:
            countGM +=bool(i & setDr)
        
       
        psym +=countCF*self.deltapsym
        psymdiff -=countCF*self.deltapsym
            
        
        
        tauDiv *=(self.capp)**countCS
        if np.any(mut[self.Tsind]==1):
            pmute=self.pmuteLOH
        tauDeath=tauDeathc+tauDiv*psymdiff
        p=(tauDeath,tauDiv*pasymdiff*(1-float(pmute)/2),tauDiv*pasymdiff*float(pmute)/2,tauDiv*psym*(1-pmute),tauDiv*psym*pmute)
        return np.asarray(p)
        
        
    
  
    
    
    
        
        
    def getPreRateMult(self,mutations,age,nc):    
         n=mutations.shape[0]
         Nc=nc*np.ones((n,1))
         P=np.zeros((n,5))
         for i in range(n):
             P[i,:]=self.getPreRate(mutations[i,:],age,Nc[i])
         return P
    
    
    
    
   
             
         
         
         




    # computes probability of group hit given mutation    
    # computes probability of group hit given mutation    
    def getProbMut(self, mutations):
        mut=mutations+np.append(1,np.zeros(17))
        
        ngrp=self.ngrp 
               
        ProbaMut=np.ones(ngrp)        
        if np.any(mut[self.Tsind]==1):
           LOHind=(mut==1).nonzero()[0]
           ProbaMut[self.Tsind]=100
           ProbaMut[LOHind]=100*1000
           ProbaMut /=ProbaMut.sum()
           
            
        else:
           ProbaMut[self.Tsind]=100/1107.0
           ProbaMut[self.Oncind]=1/1107.0
        return ProbaMut
        #increase group 1 probability if one allele hit
        
        
        
    def GetProbMutMult(self,mutations):
        ngrp=self.ngrp
        n=mutations.shape[0]
        ProbaMutMult=np.ones((n,ngrp))
        for i in range(n):
            ProbaMutMult[i,:]=self.getProbMut(mutations[i,:])
        return ProbaMutMult
    
    
    def CancerTest(self,mutations):
        
        Category1=np.logical_or(np.any(mutations[:,[2,3]]>=2,axis=1),np.any(mutations[:,[0,1]]>0,axis=1))
        Category2=np.logical_or(np.any(mutations[:,[10]]>=2,axis=1),np.any(mutations[:,np.arange(4,10)]>0,axis=1))
        
        Tscounts=np.greater_equal(mutations[:,self.Tsind],2).sum(axis=1)
        Onccounts=np.greater(mutations[:,self.Oncind],0).sum(axis=1)
        HitCounts=(Tscounts+Onccounts)>=3
        CloneSize=(self.nclones>(10**4))
        return np.any(np.logical_and(np.logical_and(HitCounts,Category1),np.logical_and(Category2,CloneSize)))
        
    
    def GetCategory(self,mutations):
        n=mutations.shape[0]
        Category=np.zeros((n,6))                               
        Category[:,0]=np.logical_or(np.any(mutations[:,[0,2,3]]>=2,axis=1),np.any(mutations[:,[1]]>0,axis=1))
        Category[:,1]=np.logical_or(np.any(mutations[:,[10]]>=2,axis=1),np.any(mutations[:,np.arange(4,10)]>0,axis=1))
        Category[:,2]=np.any(mutations[:,np.arange(11,18)]>=2,axis=1)
        Category[:,3]=np.logical_and(np.any(mutations[:,[0,2,3]]==1,axis=1),np.all(mutations[:,[1]]==0,axis=1))
        Category[:,4]=np.logical_and(np.any(mutations[:,[10]]==1,axis=1),np.all(mutations[:,np.arange(4,10)]==0,axis=1))
        Category[:,5]=np.any(mutations[:,np.arange(11,18)]==1,axis=1)
        Category=Category.astype(int)
        Category[:,[0,1,2]] *=2
        Category=np.maximum(Category[:,[0,1,2]],Category[:,[3,4,5]]).astype(int)
        v=np.dot(Category,np.array([1,3,9]).reshape((3,1)))       
        return Category,v.reshape(n)
        
    def GetChangeCategory(self,mutations,changes):
        s=changes.shape[0]
        Cat=self.GetCategory(mutations)[1]
        NewCat=self.GetCategory(np.tile(mutations,(s,1))+changes)[1]
        Ls=np.where(NewCat<>Cat)
        NewCat=NewCat[Ls]
        return NewCat,Ls[0]
    
    def Generate_Mutations(self):
        mutations=np.zeros((27,18))
        model=np.zeros((27,3))
        p=0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    
                    model[p,:]=[i,j,k]
                    p +=1
        mutations[:,[0,10,11]]=model   
        return mutations  
    
    def Generate_Mutations_Prob(self):
        Mutations_prob=np.zeros((27,27))
        G=self.Gen_mutations
        Mutations_prob_pre=self.GetProbMutMult(G)
        for i in range(27):
            C=self.GetChangeCategory(G[i,:].reshape(1,18),np.identity(18))
            for j in range(27):
                   Mutations_prob[i,j] +=Mutations_prob_pre[i,C[1][C[0]==j]].sum()
        return Mutations_prob
    
    def GeneratePreRate(self,age,nc):
        Mp=self.mutations_prob
        P=self.getPreRateMult(self.Gen_mutations,age,nc)                
        Mutations_0=np.dot(np.diag(P[:,2]),Mp)
        Mutations_1=np.dot(np.diag(P[:,4]),Mp)
        Mutations_1 +=np.diag(P[:,3])
        Death=P[:,0]
        return Mutations_0,Mutations_1,Death
    
    
    def CountCat(self,mutations,age,nclones):
        Count=np.zeros(27)
        C=self.GetCategory(mutations)[1]
        for k in range(27):
            Count[k]=nclones[np.where(C==k)[0]].sum()
        return Count.astype(int)
        
    
    def GetRate(self,age,Crypt):        
        mutations=self.mutations[self.crypt==Crypt,:]
        nclones=self.nclones[self.crypt==Crypt]
        nc=nclones.sum()
        C=self.CountCat(mutations,age,nclones)
        R=self.GeneratePreRate(age,nc)
        R_new_0=np.dot(np.diag(C),R[0])
        R_new_1=np.dot(np.diag(C),R[1])
        R_d=C*R[2]       
        return R_new_0,R_new_1,R_d,C
    
    def Get_Time_Crypt(self,age,Crypt,epsilon):
        nclones=self.nclones[self.crypt==Crypt]
        nc=nclones.sum()
        R=self.GetRate(age,Crypt)
        R_1=self.GeneratePreRate(age,nc)
        a=R[0].sum(axis=0)-R[0].sum(axis=1)-R[2]+R[1].sum(axis=0)
        b=(R[0]**2).sum(axis=0)+(R[0]**2).sum(axis=1)+R[2]**2+(R[1]**2).sum(axis=0)
        T_mu_0=np.dot(np.diag(a),R_1[0])
        T_mu_1=np.dot(np.diag(a),R_1[1])
        T_mu_d=a*R_1[2]
        T_sig_0=np.dot(np.diag(b),R_1[0])
        T_sig_1=np.dot(np.diag(b),R_1[1])
        T_sig_d=b*R_1[2]
        Rate_tot=R[0].sum()+R[1].sum()+R[2].sum()
        T_mu_0=epsilon*Rate_tot*1.0/np.abs(T_mu_0)
        T_mu_1=epsilon*Rate_tot*1.0/np.abs(T_mu_1)
        T_mu_d=epsilon*Rate_tot*1.0/np.abs(T_mu_d)
        T_sig_0=((epsilon*Rate_tot)**2)*1.0/T_sig_0
        T_sig_1=((epsilon*Rate_tot)**2)*1.0/T_sig_1
        T_sig_d=((epsilon*Rate_tot)**2)*1.0/T_sig_d
        Time=min(np.min(T_mu_0),np.min(T_mu_1),np.min(T_mu_d),np.min(T_sig_0),np.min(T_sig_1),np.min(T_sig_d))
        return Time
    
    def Get_Time_Tissue(self,age,epsilon):
        C=np.unique(self.crypt)
        T=np.zeros(C.shape[0])
        for i in range(C.shape[0]):
            T[i]=self.Get_Time_Crypt(age,C[i],epsilon)
        return T,T.min(),T.mean(),C
    
    def Simulate_Crypt(self,age):
        N=self.nc
        Cat=self.GetCategory(self.mutations)[1]
        Reactions=np.zeros((N,2))
        
        Cr=np.unique(self.crypt)
        
        for k in range(Cr.shape[0]):
            
            N_sub=np.where(self.crypt==Cr[k])[0]
            
            A=self.GetRate(age,Cr[k])[0]
            B=self.GetRate(age,Cr[k])[1]
            C=self.GetRate(age,Cr[k])[2]
            P_A=np.random.poisson(A)
            indA=np.unique(np.where((P_A>0)[0]))
            
            for i in indA:
                typ=np.where(Cat==i)[0]
                cl=np.intersect1d(N_sub,typ)
                P=(self.nclones[cl])/float(self.nclones[cl].sum())
                Death_Mut=np.random.multinomial(P_A[i,:].sum(),P)
                Reactions[cl,0] -=Death_Mut
                Reactions[cl,1] +=Death_Mut
            P_B=np.random.poisson(B)
            
            
            D_B=np.diag(P_B)
            
            
            P_B =P_B-D_B
            indB=np.unique(np.where((P_B>0)[0]))
            for i in indB:
                typ=np.where(Cat==i)[0]
                cl=np.intersect1d(N_sub,typ)
                P=(self.nclones[cl])/float(self.nclones[cl].sum())
                Birth_Mut=np.random.multinomial(P_B[i,:].sum(),P)                
                Reactions[cl,1] +=Birth_Mut 
            
               
            ind_D_B=np.where(D_B>0)[0]
            
            for i in ind_D_B:
                typ=np.where(Cat==i)[0]
                cl=np.intersect1d(N_sub,typ)
                P=(self.nclones[cl])/float(self.nclones[cl].sum())
                Birth=np.random.multinomial(D_B[i],P)                
                Reactions[cl,0] +=Birth 
            
            ind_D_B=np.where(D_B>0)[0]
             
            P_C=np.random.poisson(C)
            
            indC=np.where(P_C>0)[0]
            
            for i in indC:
                typ=np.where(Cat==i)[0]
                cl=np.intersect1d(N_sub,typ)
                P=(self.nclones[cl])/float(self.nclones[cl].sum())
                Death=np.random.multinomial(P_C[i],P)                
                Reactions[cl,0] -=Death 
        return Reactions
            
                
            
            
        
            
        
        
        
        
    
    
        
        
        
        
        
        
        
        
    
    
                   
               
    
    
        
        
   
    
    
    
            
            
            
        
        

        

        
        
        
        
        
        
        
        
    
    
        
                        
                    
                    
                    
                    

                
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

    # simulates evolution while tracking clonal expansions

    def testCancerClonal(self, silent=True):
        count=0
        t = 0
        zgrp = np.zeros(self.ngrp+1)
        for k in range(self.ngrp):
            zgrp[k+1] = zgrp[k] + self.pgrp[k]
        idMax = 1
        cancer = False
        self.mutations = np.zeros([0, self.ngrp], dtype=long)
        self.type = np.zeros(0, dtype=long)
        self.nclones = np.zeros(0, dtype=long)
        self.id = np.zeros(0, dtype=long)
        self.parent = np.zeros(0, dtype=long)
        self.birth = np.zeros(0, dtype=long)
        self.npergroup = np.zeros(self.ngrp, dtype=long)
        self.crypt = np.zeros(0, dtype=long)
        self.nc = 0
        self.nmut = 0
        self.history = np.zeros(0, dtype=object)
        #N_norm=10
        
        
        self.Times=[0]
       
        while t < self.T and not cancer:
            
            
            #N_old=self.N
            nb_crypt_old=self.nb_crypt
            
            if t<72:
                self.nb_crypt=round((1+1388*t)/100.0+1)
            elif (t>=72) and (t<20*52*2):
                self.nb_crypt=round((1+1388*72+463*(t-72))/100.0+1)
            else:
                self.nb_crypt=round((1+1388*72+463*(20*52*2-72))/100.0+1)
                
            nb_crypt_diff=self.nb_crypt-nb_crypt_old
            

            self.N=10*nb_crypt_old+10*nb_crypt_diff
          
            
            
            
            
            
            
            
            #print t
            
            
            for k in range(1, self.ngrp+1):
                self.groupsize[k] = self.nclones[self.type==k].sum()
            if self.nmut > 0:
                if not silent:
                    st = '{0:f} {1:d} ({2:d}) '.format(t/52, self.nmut, self.nclones.shape[0])
                    for k in range(0, self.ngrp):
                        st+= '-- {0:d} ({1:d}+{2:d}) '.format(self.nclones[self.mutations[:,k]>0].sum(), (self.mutations[:,k]==1).sum(),np.logical_and(self.mutations[:,k]>1,self.mutations[:,0]>1).sum())
                    st += ' -- {0:d} {1:d}'.format(sum(self.nclones > self.t3Threshold), self.nclones.max())
                    print st
                if self.CancerTest(self.mutations) :
                    cancer = True
                    break
            else:
                if not silent:
                    print t, self.nmut

            #changing
            newm0 = np.zeros([0, self.ngrp], dtype=long)
            newt0 = np.zeros(0, dtype=long)
            newc0 = np.zeros(0, dtype=long)
            newb0 = np.zeros(0, dtype=long)
            newcr0 = np.zeros(0, dtype=long)
            newpar0 = np.zeros(0, dtype=long)
            newid0 = np.zeros(0, dtype=long)
            nm = np.zeros(self.nc, dtype=long)
            idm = np.zeros(self.nc, dtype=long)
            newhist0 = np.zeros(0, dtype=object)
            #V=self.Simulate_Crypt(1)
            #print V
            for j in range(self.nc):
                nb_in_crypt = sum(self.nclones[self.crypt==self.crypt[j]])
                Poisson=np.random.poisson(self.getPreRate(self.mutations[j],1,nb_in_crypt)*0.5*self.nclones[j])
                
                self.nclones[j] +=-Poisson[0]+Poisson[3]-Poisson[2]
                if self.nclones[j]==0:
                    count=count+1
                nm[j] =Poisson[2]+Poisson[4]
                self.dcell +=Poisson[1]+Poisson[2]#######################################################################################
                idm[j] = self.id[j]
            nn = (nm).sum()
            newm = np.zeros([nn, self.ngrp], dtype=long)
            newt = np.zeros(nn, dtype=long)
            newc = np.zeros(nn, dtype=long)
            newcr = np.zeros(nn, dtype=long)
            newb = np.zeros(nn, dtype=long)
            newid = np.zeros(nn, dtype=long)
            newpar = np.zeros(nn, dtype=long)
            newhist = np.zeros(nn, dtype=object)
            ll = 0
            for j in range(self.nc):
                for k in range(nm[j]):
                    newm[ll,:] = self.mutations[j]
                    v = np.nonzero(np.random.multinomial(1,self.getProbMut(self.mutations[j,:])) > 0 )[0][0]
                    newm[ll,v] += 1
                    newt[ll] = (newm[ll,:]>0).sum()
                    newc[ll] = 1
                    newcr[ll] = self.crypt[j]
                    newb[ll] = t
                    newid[ll] = idMax
                    newpar[ll] = idm[j]
                    newhist[ll] = self.history[j] + ';' + str(t) + ' ' + str(v+1)
                    idMax += 1
                    ll += 1
            newm0 = np.concatenate((newm0, newm))
            newt0 = np.concatenate((newt0, newt), axis=0)
            newc0 = np.concatenate((newc0, newc))
            newcr0 = np.concatenate((newcr0, newcr))
            newb0 = np.concatenate((newb0, newb))
            newid0 = np.concatenate((newid0, newid))
            newpar0 = np.concatenate((newpar0, newpar))
            newhist0 = np.concatenate((newhist0, newhist))

            IS = self.nclones > 0
            self.mutations = np.concatenate((self.mutations[IS,:], newm0))
             
            
            self.type = np.concatenate((self.type[IS], newt0))
            self.nclones = np.concatenate((self.nclones[IS], newc0))
            self.crypt = np.concatenate((self.crypt[IS], newcr0))
            self.birth = np.concatenate((self.birth[IS], newb0))
            self.parent = np.concatenate((self.parent[IS], newpar0))
            self.id = np.concatenate((self.id[IS], newid0))
            self.history = np.concatenate((self.history[IS], newhist0))
            
           
           # Mutations from non exsisting clones
            
           
            R_0=self.getPreRate(np.zeros(18),1,0)
            selfcryptinter=self.crypt
            Mut_crypt=np.unique(selfcryptinter)           
            Mut_nb_cr=Mut_crypt.shape[0]
           
                
            
           #Mutations of the "normal" cells in "normal" crypts
            
            
            n = np.random.poisson(5*(nb_crypt_old-Mut_nb_cr)*R_0[np.array([2,4])].sum(),1)[0]
            newm = np.zeros([n, self.ngrp], dtype = long)
            newhist = np.zeros(n, dtype=object)
            
            
            
            
            
            for k in range(n):
                v = np.nonzero(np.random.multinomial(1,self.getProbMut(np.zeros(18))) > 0 )[0][0]
                newm[k,v] = 1
                newhist[k] = str(t) + ' ' + str(v+1)
                
                
                
                
                
            newt = np.ones(n, dtype=long)
            newc = np.ones(n, dtype=long)
            newcr =np.random.random_integers(1, nb_crypt_old, n)
            
            while (np.setdiff1d(newcr,Mut_crypt).shape[0]<n):
                   newcr =np.random.random_integers(1, nb_crypt_old, n)
                   #print newcr
            
            newpar = np.zeros(n, dtype=long)
            newb = t*np.ones(n, dtype=long)
            newid = range(idMax, idMax+n)
            idMax += n
            
            self.type = np.append(self.type, newt)
            self.parent = np.append(self.parent, newpar)
            self.birth = np.append(self.birth, newb)
            self.id = np.append(self.id, newid)
            self.mutations=np.append(self.mutations, newm, axis=0)
            
            
            self.nclones = np.append(self.nclones, newc)
            self.crypt = np.append(self.crypt, newcr)
            self.history = np.append(self.history, newhist)
            
            #print newcr
            #print newc.shape[0]
            
            
            #Mutations of nomral cells in crypts containing  mutations
            
            n_array=np.zeros(Mut_nb_cr,dtype=int)
            
            for k in range(Mut_nb_cr):
                
                size_crypt=sum(self.nclones[self.crypt==Mut_crypt[k]])
                R_0=self.getPreRate(np.zeros(18),1,size_crypt)
                n_array[k]=np.random.poisson(10*R_0[np.array([2,4])].sum(),1)[0]
                
            updated_n_array=  np.nonzero(n_array)[0]
            
            if updated_n_array.shape[0]>0:
               
               Mut_crypt_updated=Mut_crypt[updated_n_array.astype(int)]
               
               for i in range(Mut_crypt_updated.shape[0]):
                    #print n_array
                    #print updated_n_array
                    n = int(n_array[updated_n_array[i]])
                    
                    newm = np.zeros([n, self.ngrp], dtype = long)
                    newhist = np.zeros(n, dtype=object)
                    for k in range(n):
                         v = np.nonzero(np.random.multinomial(1,self.getProbMut(np.zeros(18))) > 0 )[0][0]
                         newm[k,v] = 1
                         newhist[k] = str(t) + ' ' + str(v+1)
                   
                    newt = np.ones(n, dtype=long)
                    newc = np.ones(n, dtype=long)
                    newcr =Mut_crypt_updated[i]*np.ones(n,dtype=int)
                    #print n
                    #print newcr
                    newpar = np.zeros(n, dtype=long)
                    newb = t*np.ones(n, dtype=long)
                    newid = range(idMax, idMax+n)
                    idMax += n
                    
                    
                    #print newcr
                    #print newc.shape[0]
            
                    self.type = np.append(self.type, newt)
                    self.parent = np.append(self.parent, newpar)
                    self.birth = np.append(self.birth, newb)
                    self.id = np.append(self.id, newid)
                    self.mutations=np.append(self.mutations, newm, axis=0)
            
            
                    self.nclones = np.append(self.nclones, newc)
                    self.crypt = np.append(self.crypt, newcr)
                    self.history = np.append(self.history, newhist)
            
                   
                    self.nc = self.nclones.shape[0]
                    
                    
                    
           #Specific to the Growth phase         
           
            if (t<52*20):
               Mut_crypt=np.unique(self.crypt)           
               Mut_nb_cr=Mut_crypt.shape[0]
               size_crypt_array=np.zeros(Mut_nb_cr)
               for i in range(Mut_nb_cr):
                   size_crypt_array[i]=sum(self.nclones[self.crypt==Mut_crypt[k]])+10
               P=np.array([size_crypt_array.sum(),(nb_crypt_old-Mut_nb_cr)*10])
               P=P/float(P.sum())
               Cr_daughters=np.random.multinomial(nb_crypt_diff,P)[0]
               if Cr_daughters>0:
                  size_crypt_array=size_crypt_array/float(size_crypt_array.sum())
                  Cr_daughters_dist=np.random.multinomial(Cr_daughters,size_crypt_array)
                  Cr_daughters_mutant=np.nonzero(Cr_daughters_dist)[0]
                  cr_count=int(nb_crypt_old)
                  #print Cr_daughters_mutant
                  
                  for i in range(Cr_daughters_mutant.shape[0]):
                      ind=Cr_daughters_mutant[i]
                      #print ind
                     
                      cr_count +=1
                      crypt=Mut_crypt[ind]
                      copies=Cr_daughters_dist[ind]
                      #print copies
                      indices=np.where(self.crypt==crypt)[0]
                   
                      newhist=np.zeros(indices.shape[0],dtype=object)
                   
                      for k in range(indices.shape[0]):
                          newhist[k]=self.history[indices[k]]+';'+str(t)+ ' ' +'split'
                       
                   
                      newm=self.mutations[indices,:]
                      newt=self.type[indices]
                      newc=self.nclones[indices]
                      newcr=cr_count*np.ones(newc.shape[0],dtype=int)
                      newpar=indices
                      newb=t*np.ones(indices.shape[0])
                      newid = range(idMax, idMax+copies*indices.shape[0])
                      idMax +=copies*indices.shape[0]
                      
                      
                      #print newcr
                      #print newc.shape[0]
                      
                      
                      self.type = np.append(self.type, npmat.repmat(newt,1,copies)[0])
                      self.parent = np.append(self.parent, npmat.repmat(newpar,1,copies)[0])
                      
                      self.birth = np.append(self.birth, npmat.repmat(newb,1,copies)[0])
                      self.id = np.append(self.id, newid)
                      self.mutations=np.append(self.mutations, npmat.repmat(newm,copies,1), axis=0)
            
                      #print npmat.repmat(newcr,1,copies)[0]
                      #print copies
                      self.nclones = np.append(self.nclones,npmat.repmat(newc,1,copies)[0])
                      self.crypt = np.append(self.crypt, npmat.repmat(newcr,1,copies)[0])
                      self.history = np.append(self.history, npmat.repmat(newhist,1,copies)[0])
            
            self.nc = self.nclones.shape[0]
                       
                       
            
            
            #print N_norm
            if t>1000 and (t % 100)==0:
               self.Times=np.append(self.Times,self.Get_Time_Tissue(1,0.1)[2])
            t = t+self.delta
            

            self.nc = self.nclones.shape[0]
            self.nmut = self.nclones.sum()
            for k in range(self.ngrp):
                self.npergroup[k] = self.nclones[self.mutations[:,k]>0].sum()            

        st = 'Time {0:f} Total size {1:d} ({2:d})   Groups '.format(t/52, self.nmut, self.nclones.shape[0])
        for k in range(0, self.ngrp):
            st+= '-- {0:d} : {1:d} ({2:d}) '.format(k+1, self.nclones[self.mutations[:,k]>0].sum(), (self.mutations[:,k]>0).sum())
        st += ' -- {0:d}'.format(sum(self.nclones > self.t3Threshold))
        print st
        st = '     Levels '
        for k in range(0, self.ngrp):
            st+= '-- {0:d}: {1:d} ({2:d}) '.format(k+1, self.nclones[self.type==k+1].sum(), (self.type==k+1).sum())
        print st
        st = '     Large Clones {0:d} '.format((self.nclones > self.t3Threshold).sum())
        print st
        self.lastTime = t
        self.cancer = cancer

        return  [cancer,(self.nclones > self.t3Threshold).sum(),self.lastTime,self.nclones[self.type==1].sum(),self.nclones[self.type==2].sum(),self.nclones[self.type==3].sum(),self.nclones[self.type==4].sum(),self.nclones[self.type==5].sum(),count]

























count=0

tab6=[]
tabcr=[]
tabtime=[]
tabn=[]
tabncancer=[]
tabhistory=[]
tabmutations=[]
a=False
Times=np.zeros((50,30))
for i in range(100):
    print i
    t=tissue()
    [a,b,c,a1,a2,a3,a4,a5,a6]=t.testCancerClonal()
    n=sum(t.nclones[t.crypt==t.crypt[np.argmax(t.nclones)]])
    tabcr.append(n)
    #Times[i,range(t.Times.shape[0])]=t.Times
    if a:
    #tab6.append(float(a6)/b)
      
      tabtime.append(c)
      tabmutations.append(t.mutations)
      tabhistory.append(t.history)
      tabncancer.append(n)
    
    tabn.append(n)
   
    if ((t.type>1).sum())>0:
        count +=1
    print n
tabcr.sort()
print "Tabnc is :"+str (tabcr)
print "Tabtime is :"+str(tabtime)
print "Tabhistory is :"+str(tabhistory)
print "Tabmutations is :"+str(tabmutations)
print "Tabncancer is :"+str(tabncancer)
print "Tabtime is:"+str(Times)
   
        
