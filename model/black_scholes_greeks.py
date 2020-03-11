
# coding: utf-8

# In[17]:


import math
from scipy.stats import norm as n

# Compute d1,d2,expire time 
def probability(S,K,r,q,t,T,sigma):
    time_value = T-t
    d1 = (math.log(S/K) + time_value*(r-q+math.pow(sigma,2)/2))/ (sigma*math.sqrt(time_value))
    d2 = d1 - sigma*math.sqrt(time_value)
    return [time_value,d1,d2]

# Compute call Option and put Option 
def callPutOptionPrices(S,K,r,q,t,T,sigma):
      # get value of d1,d2 and expiry value

    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    time_value,d1,d2 = probability(S,K,r,q,t,T,sigma) 

    # time_value,d1,d2 = probability_list[0],probability_list[1],probability_list[2]
    # d1 = probability_list[1]
    # d2 = probability_list[2]
    c = (S*math.exp(-q*time_value)*n.cdf(d1)) - (K*math.exp(-r*time_value)*n.cdf(d2))
    p = (K*math.exp(-r*time_value)*n.cdf(-d2)) - (S*math.exp(-q*time_value)*n.cdf(-d1))
    return [c,p]

# Compute delta
def deltas(S,K,r,q,t,T,sigma):
     # get value of d1,d2 and expiry value

    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    time_value,d1,d2 = probability(S,K,r,q,t,T,sigma) 


    # time_value,d1 = probability_list[0],probability_list[1]
    # only retrieve d1 as d2 is not needed for delta calculation
    # d1 = probability_list[1]
    #d2 = probability_list[2]
    
    deltaC = math.exp(-q*time_value)*n.cdf(d1)
    deltaP = -math.exp(-q*time_value)*n.cdf(-d1)
    
    deltaC100 = 100*deltaC
    deltaP100 = 100*deltaP
    
    return [deltaC, deltaP, deltaC100, deltaP100]

#compute lambdas

## Aaron - If lambdas function is not used in the black scholes function below,
##         can we just get rid of it?
##         However there is a lesson point that is included in lambdas.

def lambdas (S,K,r,q,t,T,sigma):

    #Lesson - "_" is used as a filler for unused variables
    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    delta_output1, delta_output2 , _ , _ = deltas(S,K,r,q,t,T,sigma)
    callPut_output1 , callPut_output2 = callPutOptionPrices(S,K,r,q,t,T,sigma)
    lambdasC = delta_output1 * S/callPut_output1
    lambdasP = delta_output2 * S/callPut_output2


    # delta_output = deltas(S,K,r,q,t,T,sigma)
    # callPut_output = callPutOptionPrices(S,K,r,q,t,T,sigma)
    # lambdasC = delta_output[0] * S/callPut_output[0]
    # lambdasP = delta_output[1] * S/callPut_output[1]
    return [lambdasC, lambdasP]


def blackScholes(S,K,r,q,t,T,sigma):

    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    time_value,d1,d2 = probability(S,K,r,q,t,T,sigma) 
    # probability_list = probability(S,K,r,q,t,T,sigma) 
    # time_value = probability_list[0]
    # d1 = probability_list[1]
    # d2 = probability_list[2]
    


    c,p = callPutOptionPrices(S,K,r,q,t,T,sigma)
#     c = callPut_output[0]
#     p = callPut_output[1]
    
    deltaC,deltaP,_,_ = deltas(S,K,r,q,t,T,sigma)
#     deltaC = delta_output[0]
#     deltaP = delta_output[1]
    
    # lambda_output = lambdas (S,K,r,q,t,T,sigma)
    # lambdasC = lambda_output[0]
    # lambdasP = lambda_output[1]
    
    # callPut_output = callPutOptionPrices(S,K,r,q,t,T,sigma)
    # c = callPut_output[0]
    # p = callPut_output[1]
    
    # delta_output = deltas(S,K,r,q,t,T,sigma)
    # deltaC = delta_output[0]
    # deltaP = delta_output[1]
    
    # lambda_output = lambdas (S,K,r,q,t,T,sigma)
    # lambdasC = lambda_output[0]
    # lambdasP = lambda_output[1]
    

    gammaC = math.exp(-q*time_value) * ((n.pdf(d1))/(S*sigma*math.sqrt(time_value)))
    gammaP = gammaC
    
    vegaC = S*math.exp(-q*time_value) * n.pdf(d1) * math.sqrt(time_value)*0.01
    vegaP = vegaC
    
    rhoC = K*time_value*math.exp(-r*time_value) * n.cdf(d2)*0.01
    rhoP = -K*time_value*math.exp(-r*time_value) * n.cdf(-d2)*0.01
    
    #Lesson - Implicit joining 

    psiC = ((-S*time_value*math.exp(-q*time_value)*n.cdf(d1) - 
        (S*math.sqrt(time_value)/sigma)*math.exp(-q*time_value)  
        * n.pdf(d1)+(K*math.sqrt(time_value)/sigma)*math.exp(-r*time_value) * n.pdf(d2)) * 0.01 )

    psiP = ((S*time_value*math.exp(-q*time_value)*n.cdf(-d1) - 
        (S*math.sqrt(time_value)/sigma)*math.exp(-q*time_value) 
        * n.pdf(-d1)+(K*math.sqrt(time_value)/sigma)*math.exp(-r*time_value) * n.pdf(-d2)) * 0.01 )
    
    SSc = (c-S*deltaC)/K
    SSp = (p-S*deltaP)/K
    
    IVc = max(S-K, 0)
    IVp = max(K-S, 0)
    
    TVc = c-IVc
    TVp = p-IVp
    
    ZVc = max(S*math.exp((r-q)*time_value)-K, 0) * math.exp(-r*time_value)
    ZVp = max(K-S*math.exp((r-q)*time_value), 0) * math.exp(-r*time_value)
    
    # delta,c,p,lambdas not return in the above function
    return (gammaC,gammaP,vegaC,vegaP,rhoC,rhoP,psiC,psiP,SSc,SSp,IVc,IVp,TVc,TVp,ZVc,ZVp)
    
# compute gamma1%, placed below blacksholes function as it's function not needed for the calculation
def gamma1Percent(S,K,r,q,t,T,sigma):



    #Lesson - "_" is used as a filler for unused variables
    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    delta_output1,delta_output2 , _ ,_ = deltas(S*1.01,K,r,q,t,T,sigma)
    delta_koutput1,delta_koutput2 , _ , _ = deltas(S*0.99,K,r,q,t,T,sigma)

    GammaC = (delta_output1-delta_koutput1)/2
    GammaP = (delta_output2-delta_koutput2)/2
    

#     delta_output1 = deltas(S*1.01,K,r,q,t,T,sigma)
#     delta_output2 = deltas(S*0.99,K,r,q,t,T,sigma)

#     GammaC = (delta_output1[0]-delta_output2[0])/2
#     GammaP = (delta_output1[1]-delta_output2[1])/2
    

    return(GammaC,GammaP)

# compute theta1Day
def theta1Day(S,K,r,q,t,T,sigma):
   

    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    Vc1,Vc2 = callPutOptionPrices(S,K,r,q,t+1/365,T,sigma)
    Vck1,Vck2 = callPutOptionPrices(S,K,r,q,t,T,sigma)
    
    thetaC1 = Vc1-Vck1
    thetaP1 = Vc2-Vck2
   

#     Vc1 = callPutOptionPrices(S,K,r,q,t+1/365,T,sigma)
#     Vc2 = callPutOptionPrices(S,K,r,q,t,T,sigma)
    
#     thetaC1 = Vc1[0]-Vc2[0]
#     thetaP1 = Vc1[1]-Vc2[1]
   
    return(thetaC1,thetaP1)

# compute theta7Day
def theta7Days(S,K,r,q,t,T,sigma):
    

    #Lesson - Multiple assignment of variables
    #Lesson - Assigning variables straight to function returns
    callPut_output1,callPut_output2 = callPutOptionPrices(S,K,r,q,t+7/365,T,sigma)
    callPut_outputk1,callPut_outputk2 = callPutOptionPrices(S,K,r,q,t,T,sigma)
    
    thetaC7 = callPut_output1-callPut_outputk1
    thetaP7 = callPut_output2-callPut_outputk2
   

    # callPut_output1 = callPutOptionPrices(S,K,r,q,t+7/365,T,sigma)
    # callPut_output2 = callPutOptionPrices(S,K,r,q,t,T,sigma)
    
    # thetaC7 = callPut_output1[0]-callPut_output2[0]
    # thetaP7 = callPut_output1[1]-callPut_output2[1]
   
    return(thetaC7,thetaP7)
    
# print(callPutOptionPrices(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4))
# print(deltas(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4)) 
# print(blackScholes(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4))
# print(gamma1Percent(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4))
# print(theta1Day(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4))
# print(theta7Days(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4))

   
    


# In[8]:


#!pip3 install numpy==1.17.0

