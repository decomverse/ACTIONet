import pandas as pd
import random 
import numpy as np 
from anndata import AnnData
from typing import Optional
import scipy


def string_list_to_int_list(string_list):
    '''
    maps strings to ints
    i.e. ['a','b','c','b','a']--> [1,2,3,2,1]
    '''
    #generate the map
    string_to_int={}
    index=0
    for entry in string_list:
        if entry not in string_to_int:
            string_to_int[entry]=index
            index+=1
    int_list=[0]*len(string_list)
    for i in range(len(string_list)):
        int_list[i]=string_to_int[string_list[i]]
    return int_list


def compute_phi(A,labels,s0,s1,s2):
    n=len(labels)
    #handle labels passed as list by wrapping them in a pandas series object
    if type(labels) != type(pd.Series): 
        labels=pd.Series(labels)
    counts=labels.value_counts()
    categories=counts.index.tolist()
    w=A.data 
    pvec=counts/n
    k=len(pvec)
    m1_rawphi=s0/(n*(n-1))*(n**2*k*(2-k)-n*sum(1/pvec))
    Q1 = sum(1/pvec)
    Q2 = sum(1/np.power(pvec,2))
    Q3 = sum(1/np.power(pvec,3))
    Q22=np.sum( np.expand_dims(1/pvec,axis=1)*np.expand_dims(1/pvec,axis=0))
    E1 = (np.power(n,2) * Q22 - n * Q3)/(n * (n - 1))
    E2 = 4 * np.power(n,3) * Q1 - 4 * np.power(n,3) * k * Q1 + np.power(n,3) * np.power(k,2) * Q1 - 2 * (2 * np.power(n,2) * Q2 - np.power(n,2) * k * Q2) + 2 * n * Q3 - np.power(n,2) * Q22
    E2 = E2/(n * (n - 1) * (n - 2))
    A1 = 4 * np.power(n,4) * np.power(k,2) - 4 * np.power(n,4) * np.power(k,3) + np.power(n,4) * np.power(k,4) - (2 * np.power(n,3) * k * Q1 - np.power(n,3) * np.power(k,2) *Q1)
    A2 = 4 * np.power(n,3) * Q1 - 4 * np.power(n,3) * k * Q1 + np.power(n,3) * np.power(k,2) * Q1 - (2 * np.power(n,2) * Q2 - np.power(n,2) *k * Q2)
    Apart = A1 - 2 * A2
    B1 = 4 * np.power(n,3) * Q1 - 4 * np.power(n,3) * k * Q1 + np.power(n,3) * np.power(k,2) * Q1 - (2 * np.power(n,2) * Q2 - np.power(n,2) * k * Q2)
    B2 = 2 * np.power(n,2) * Q2 - np.power(n,2) * k * Q2 - n * Q3
    B3 = np.power(n,2) * Q22 - n * Q3
    Bpart = B1 - B2 - B3
    C1 = 2 * np.power(n,3) * k * Q1 - np.power(n,3) * np.power(k,2) * Q1 - np.power(n,2) * Q22
    C2 = 2 * np.power(n,2) * Q2 - np.power(n,2) * k * Q2 - n * Q3
    Cpart = C1 - 2 * C2
    E3 = (Apart - 2 * Bpart - Cpart)/(n * (n - 1) * (n - 2) * (n - 3))
    m2_rawphi = s1 * E1 + (s2 - 2 * s1) * E2 + (np.power(s0,2) - s2 + s1) * E3
    v_i=labels[A.row]
    v_j=labels[A.col]
    p_i = np.asarray(pvec[v_i])
    p_j = np.asarray(pvec[v_j])
    rawphi = int(sum(w * (2 * (v_i.reset_index(drop=True) == v_j.reset_index(drop=True)) - 1)/(p_i * p_j)))
    mean_rawphi = m1_rawphi
    var_rawphi = m2_rawphi - np.power(mean_rawphi,2)
    phi_z = (rawphi - mean_rawphi)/np.sqrt(var_rawphi)
    phi_logPval=-1*np.log10(scipy.stats.norm.sf(phi_z))
    dictz=phi_z
    logPval=phi_logPval
    phi=rawphi
    return dictz, logPval, phi


def assess_categorical_autocorrelation(
        adata: AnnData,
        labels: list,
        perm_no: int=100):
    A=adata.obsp['ACTIONet']
    #if labels is a list of strings, change them to a numerical encoding
    classnames, labels = np.unique(labels, return_inverse=True)
    labels=labels+1
    
    #labels=string_list_to_int_list(labels)
    w=A.data
    s0=sum(w)
    s1=sum(4*w**2)/2
    s2=int(sum(np.power(A.sum(axis=1)+A.sum(axis=0).transpose(),2)))
    A=A.tocoo()
    dictz,logPval,phi=compute_phi(A,labels,s0,s1,s2)
    rand_phis=[]
    #set the random seed
    np.random.seed(0)
    random.seed(0)    
    for i in range(perm_no):
        rand_labels=random.sample(labels.tolist(),labels.shape[0])
        rand_dictz, rand_logPval, rand_phi=compute_phi(A,rand_labels,s0,s1,s2)
        rand_phis.append(rand_phi)
    z=(phi-np.mean(rand_phis))/np.std(rand_phis)
    return z

