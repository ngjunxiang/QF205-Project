import numpy as np


import black_scholes_lib

# Lesson - importing functions from other files
# Lesson - assigning a function to a variable name
checker = black_scholes_lib.check_params

def blackScholes_cranknicolson(S,K,r,q,t,T,sigma,M,N):

    #Lesson - variable assigned function
    if not checker(S, K, r, q, T, sigma, M, N):
        return None

    #initializing values Smax, deltaT and deltaS based on assumptions (S<2K)
    Smax = 2*K
    deltaT = T/N
    deltaS = Smax/M
    
    # Initialize fCallOption and fPutOption
    fCallOption = np.zeros((N+1,M+1))
    fPutOption = np.zeros((N+1,M+1))
    
    for j in range(M+1):
        fCallOption[N,j] = max((j*deltaS)-K,0)
        fPutOption[N,j] = max((K-(j*deltaS)),0)
    
    
    # Creating constants alpha j, beta j and gamma j
    
#     def alphaj(j):
#         aj = 0.25 * deltaT * (sigma**2 * j**2 - (r-q) * j)
#         return aj 
    # Lesson - Lambda expressions
    alphaj = lambda x: 0.25 * deltaT * (sigma**2 * j**2 - (r-q) * j)
    
    
#     def betaj(j):
#         bj = -0.5 * deltaT * (sigma**2 * j**2 + r)
#         return bj 
    
    # Lesson - Lambda expressions
    betaj = lambda x: -0.5 * deltaT * (sigma**2 * j**2 + r)

#     def gammaj(j):
#         gj = 0.25 * deltaT * (sigma**2 * j**2 + (r-q) * j)
#         return gj
    # Lesson - Lambda expressions
    gammaj = lambda x: 0.25 * deltaT * (sigma**2 * j**2 + (r-q) * j)


    # Initialise matrix M2
    M2 = np.zeros((M+1,M+1))
    M2[0,0],M2[M,M] = 1,1
    
    #Initialise Matrix M1
    M1=np.zeros((M+1,M+1))
    M1[0,0],M1[M,M] = 1,1
    
    # Populate matrix M1 and M3 with alpha j, beta j, gamma j
    for j in range(1,(M-1)+1):
        
        # Assignment of values to positions with regards to alpha j
        M2[j,j-1] = alphaj(j) 
        M1[j,j-1] = -1 * alphaj(j)
        
        # Assignment of values to positions with regards to beta j
        M2[j,j] = 1+betaj(j)
        M1[j,j] = 1-betaj(j)
        
        # Assignment of values to positions with regards to gamma j
        M2[j,j+1] =gammaj(j)
        M1[j,j+1] =-1* gammaj(j)
    
    
    inverse_matrixM1 = np.linalg.inv(M1)
    
    for i in range(N-1, 0-1, -1):
        # Transpose every row in fCallOption from horizontal to vertical for matrix multiplication with M2
        Fi_addone_call_option = np.transpose(fCallOption[[i+1]])
        Fi_addone_put_option = np.transpose(fPutOption[[i+1]])
        
        # Matrix multiplication of M2 and call option
        b_hat_call = M2 @ Fi_addone_call_option
        b_hat_put = M2 @ Fi_addone_put_option
         
        
        # Initialize bCall and bPut
        b_call = np.zeros((M+1,1)) #M+1 because j=0,1...M so M+1
        b_put = np.zeros((M+1,1))
        
		#Populate bCall and bPut
        # Lesson - Explicit joining
        b_call[0,0] = 0
        b_call[M,0] = Smax*np.exp((-1*q)*(N-i)*(deltaT)) -   \
                        K*np.exp((-1*r)*(N-i)*(deltaT)) # not sure if this is right should be Smax - K*np.exp((-1*r)*(N-i)*(deltaT) according to the formula sheet
        
        b_put[0,0] = K*np.exp((-1*r)*(N-i)*(deltaT))
        b_put[M,0] = 0
        
        for j in range(1,(M-1)+1):
            b_call[j] = b_hat_call[j]
            b_put[j] = b_hat_put[j]
        
		# Calculate Fi for call option and put option
        Fi_call_option = inverse_matrixM1 @ b_call
        Fi_put_option = inverse_matrixM1 @ b_put
	
        
        # Utilize previous fCallOption and fPutOption containers to store new Fi values.
        for j in range(0, M+1):
            fCallOption[i,j] = Fi_call_option[j,0]
            fPutOption[i,j] = Fi_put_option[j,0]
			
        
    # Calculate k
    k = int(np.floor(S/deltaS))
    
    # Calculate Call option prices and put option prices

    # Lesson - Explicit joining
    call_option_price = fCallOption[0,k] + ((fCallOption[0,k+1] - \
        fCallOption[0,k])/deltaS)*(S - (k * deltaS))

    put_option_price = fPutOption[0,k] + ((fPutOption[0,k+1] - \
        fPutOption[0,k])/deltaS)*(S - (k * deltaS))
    
    return call_option_price, put_option_price

# print(blackScholes_cranknicolson(S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183/365), sigma=0.4, M=100, N=10_000))
    