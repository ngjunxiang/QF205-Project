
# coding: utf-8

# In[17]:


import math
from scipy.stats import norm as n

# Compute d1,d2,expire time 
def probability(S,K,r,q,T,sigma,t):
    time_value = T-t
    d1 = (math.log(S/K) + time_value*(r-q+math.pow(sigma,2)/2))/ (sigma*math.sqrt(time_value))
    d2 = d1 - sigma*math.sqrt(time_value)
    return [time_value,d1,d2]

# Compute call Option and put Option 
def callPutOptionPrices(S,K,r,q,T,sigma,t=0):

    # get value of d1,d2 and expiry value
    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    time_value,d1,d2 = probability(S,K,r,q,T,sigma,t) 

    c = (S*math.exp(-q*time_value)*n.cdf(d1)) - (K*math.exp(-r*time_value)*n.cdf(d2))
    p = (K*math.exp(-r*time_value)*n.cdf(-d2)) - (S*math.exp(-q*time_value)*n.cdf(-d1))
    return [c,p]

# Compute delta
def deltas(S,K,r,q,T,sigma,t=0):

    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    time_value,d1,_ = probability(S,K,r,q,T,sigma,t) 
    
    deltaC = math.exp(-q*time_value)*n.cdf(d1)
    deltaP = -math.exp(-q*time_value)*n.cdf(-d1)
    
    deltaC100 = 100*deltaC
    deltaP100 = 100*deltaP
    
    return [deltaC, deltaP, deltaC100, deltaP100]

#compute lambdas
def lambdas (S,K,r,q,T,sigma,t=0):

    #Lesson = Assigning tuple to function calling with *Args
    args = (S,K,r,q,T,sigma,t)

    #Lesson - "_" is used as a filler for unused variables
    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    #Lesson - Assigning tuple to function calling with *Args
    delta_output1, delta_output2 , _ , _ = deltas(*args)
    callPut_output1 , callPut_output2 = callPutOptionPrices(*args)
    #compute lambdas
    lambdasC = delta_output1 * S/callPut_output1
    lambdasP = delta_output2 * S/callPut_output2

    return [lambdasC, lambdasP]


def blackScholes(S,K,r,q,T,sigma,t=0):

    #Lesson = Assigning tuple to function calling with *Args
    args = (S,K,r,q,T,sigma,t)

    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    #Lesson = Assigning tuple to function calling with *Args
    time_value,d1,d2 = probability(*args) 

    #Retrieve CallPut
    c,p = callPutOptionPrices(*args)
    
    #Retrieve Delta
    deltaC,deltaP,_,_ = deltas(*args)
    
    #compute gamma
    gammaC = math.exp(-q*time_value) * ((n.pdf(d1))/(S*sigma*math.sqrt(time_value)))
    gammaP = gammaC
    
    #compute vega
    vegaC = S*math.exp(-q*time_value) * n.pdf(d1) * math.sqrt(time_value)*0.01
    vegaP = vegaC
    
    #Compute rho
    rhoC = K*time_value*math.exp(-r*time_value) * n.cdf(d2)*0.01
    rhoP = -K*time_value*math.exp(-r*time_value) * n.cdf(-d2)*0.01
    
    #Lesson - Implicit joining 
    #psi
    psiC = ((-S*time_value*math.exp(-q*time_value)*n.cdf(d1) - 
        (S*math.sqrt(time_value)/sigma)*math.exp(-q*time_value)  
        * n.pdf(d1)+(K*math.sqrt(time_value)/sigma)*math.exp(-r*time_value) * n.pdf(d2)) * 0.01 )

    psiP = ((S*time_value*math.exp(-q*time_value)*n.cdf(-d1) - 
        (S*math.sqrt(time_value)/sigma)*math.exp(-q*time_value) 
        * n.pdf(-d1)+(K*math.sqrt(time_value)/sigma)*math.exp(-r*time_value) * n.pdf(-d2)) * 0.01 )
    
    #Strike sensitive
    SSc = (c-S*deltaC)/K
    SSp = (p-S*deltaP)/K
    
    #Intrinsic value
    IVc = max(S-K, 0)
    IVp = max(K-S, 0)
    
    #Time value
    TVc = c-IVc
    TVp = p-IVp
    
    #Zero volatility
    ZVc = max(S*math.exp((r-q)*time_value)-K, 0) * math.exp(-r*time_value)
    ZVp = max(K-S*math.exp((r-q)*time_value), 0) * math.exp(-r*time_value)
    
    return (gammaC,gammaP,vegaC,vegaP,rhoC,rhoP,psiC,psiP,SSc,SSp,IVc,IVp,TVc,TVp,ZVc,ZVp)
    
# compute gamma1%, placed below blacksholes function as it's function not needed for the calculation
def gamma1Percent(S,K,r,q,T,sigma,t=0):

    #Lesson - "_" is used as a filler for unused variables
    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    delta_output1,delta_output2 , _ ,_ = deltas(S*1.01,K,r,q,T,sigma,t)
    delta_koutput1,delta_koutput2 , _ , _ = deltas(S*0.99,K,r,q,T,sigma,t)

    GammaC = (delta_output1-delta_koutput1)/2
    GammaP = (delta_output2-delta_koutput2)/2
    
    return(GammaC,GammaP)

# compute theta1Day
def theta1Day(S,K,r,q,T,sigma,t=0):
   
    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    Vc1,Vc2 = callPutOptionPrices(S,K,r,q,T,sigma,t+1/365)
    Vck1,Vck2 = callPutOptionPrices(S,K,r,q,T,sigma,t)
    
    thetaC1 = Vc1-Vck1
    thetaP1 = Vc2-Vck2
   
    return(thetaC1,thetaP1)

# compute theta7Day
def theta7Days(S,K,r,q,T,sigma,t=0):
    
    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    callPut_output1,callPut_output2 = callPutOptionPrices(S,K,r,q,T,sigma,t+7/365)
    callPut_outputk1,callPut_outputk2 = callPutOptionPrices(S,K,r,q,T,sigma,t)
    
    thetaC7 = callPut_output1-callPut_outputk1
    thetaP7 = callPut_output2-callPut_outputk2
   
    return(thetaC7,thetaP7)
    
# print(callPutOptionPrices(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4))
# print(deltas(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4)) 
# print(blackScholes(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4))
# print(gamma1Percent(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4))
# print(theta1Day(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4))
# print(theta7Days(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4))

   
    


# In[8]:


#!pip3 install numpy==1.17.0

