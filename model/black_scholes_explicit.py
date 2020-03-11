import math
import numpy as np

# Input Parameters
# S=Stock Price, K=Exercise Price, r=Interest Rate(%), q=Yiled Rate(%), t=0, T=Time to Maturity, sigma=Volatility(%), M=100, N=10_000

import black_scholes_lib

# Lesson - importing functions from other files
# Lesson - assigning a function to a variable name
# Lesson - Global variable
checker = black_scholes_lib.check_params

def black_scholes_explicit(S, K, r, q, T, sigma, M, N, t=0):

    # Check assumption
    #Lesson - variable assigned function
    if not checker(S, K, r, q, T, sigma, M, N):
        return None

    # 1) Compute dt, ds
    S_max = 2 * K
    dt = T / N
    ds = S_max / M

    # 2) Compute (Call Option) fc & (Put Option) fp

    # Setup Matrix
    fc = np.zeros((N + 1, M + 1))
    fp = np.zeros((N + 1, M + 1))

    # Compute and insert base fc = max(j(ds) - K, 0), fp = max(K - j(ds),0), for j = 0,1 ..., M
    for j in range(M + 1):
        fc[N, j] = max((j * ds) - K, 0)
        fp[N, j] = max((K - (j * ds)), 0)

    # 3) Compute Vector Fhci, Fhpi, Fic, Fip

    # Initialising matrix A
    matrix_A = np.zeros((M + 1, M + 1))
    # Inserting inital values for Matrix A
    matrix_A[0, 0], matrix_A[M, M] = 1, 1

    # Looping and inserting the values as given by algorithm
    # aj
    # bj
    # cj
    for j in range(1, M - 1 + 1, 1):
        matrix_A[j, j - 1] = 1 / 2 * dt * (sigma ** 2 * j ** 2 - (r - q) * j)
        matrix_A[j, j] = 1 - dt * (sigma ** 2 * j ** 2 + r)
        matrix_A[j, j + 1] = 1 / 2 * dt * (sigma ** 2 * j ** 2 + (r - q) * j)

    # 3.1) Compute Fhci, Fhpi
    for i in range(N - 1, 0 - 1, -1):
        Fhci = matrix_A @ fc[i + 1]
        Fhpi = matrix_A @ fp[i + 1]

        # 3.2) Compute Fic, Fip
        # Initialise Fic
        Fic = np.zeros(M + 1)
        Fic[M] = S_max - (K * math.exp(-r * (N - i) * dt))

        # Initialise Fip
        Fip = np.zeros(M + 1)
        Fip[0] = K * math.exp(-r * (N - i) * dt)

        # Adding both Fhci and Fhpi values into matrix
        for j in range(1, M - 1 + 1):
            Fic[j] = Fhci[j]
            Fip[j] = Fhpi[j]

        # Adding the computed price values into fc and fp
        for h in range(0, M + 1):
            fc[i, h] = Fic[h]
            fp[i, h] = Fip[h]

        # 4) Find K
        k = int(np.floor(S / ds))

        # 5) Option Price
        Vc = fc[0, k] + (fc[0, k + 1] - fc[0, k] / ds) * (S - K * ds)
        Vp = fp[0, k] + (fp[0, k + 1] - fp[0, k] / ds) * (S - K * ds)

    # print("Call:", Vc, "Put:", Vp)

    return (Vc, Vp)



# print(
#     black_scholes_explicit(
#         S=50.0, K=50.0, r=0.04, q=0.01, t=0, T=(183 / 365), sigma=0.4, M=100, N=10_000
#     )
# )

