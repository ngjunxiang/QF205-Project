import numpy as np

import black_scholes_lib

# Lesson - importing functions from other files
# Lesson - assigning a function to a variable name
checker = black_scholes_lib.check_params

def black_scholes_implicit(S, K, r, q, T, sigma, M, N, t=0):
    # Check assumption

    #Lesson - variable assigned function
    if not checker(S, K, r, q, T, sigma, M, N):
        return None

    # Defining Values
    S_max = 2 * K
    dt = T / N
    ds = S_max / M

    # Initialising put and call matrix
    fc = np.zeros((N + 1, M + 1))
    fp = np.zeros((N + 1, M + 1))

    # Inserting the default values for the put and call matrix
    for i in range(M + 1):
        fc[N, i] = max((i * ds) - K, 0)
        fp[N, i] = max((K - (i * ds)), 0)

    # Initialising matrix A
    matrix_A = np.zeros((M + 1, M + 1))
    # Inserting inital values for Matrix A
    matrix_A[0, 0], matrix_A[M, M] = 1, 1

    # Looping and inserting the values as defined by the algorithm
    # aj
    # bj
    # cj
    for i in range(1, M - 1 + 1, 1):
        matrix_A[i, i - 1] = 1 / 2 * dt * ((r - q) * i - sigma ** 2 * i ** 2)
        matrix_A[i, i] = 1 + dt * (sigma ** 2 * i ** 2 + r)
        matrix_A[i, i + 1] = -1 / 2 * dt * (sigma ** 2 * i ** 2 + (r - q) * i)

    # Inverse Matrix A to get A^{-1}
    inverse_A = np.linalg.inv(matrix_A)

    # Loop to Compute Fi
    for i in range(N - 1, 0 - 1, -1):
        # Initialising both f hat (call & put)
        # call
        Fhc = np.zeros((M + 1, 1))
        Fhc[0] = 0
        Fhc[M] = (S_max - (K * np.exp(-r * (N - i) * dt)))
        # put
        Fhp = np.zeros((M + 1, 1))
        Fhp[0] = (K * np.exp(-r * (N - i) * dt))
        Fhp[M] = 0

        # Looping to get the middle values for f hats
        for j in range(1, M - 1 + 1, 1):
            # pump all the values to Fh
            Fhc[j] = fc[i + 1, j]
            Fhp[j] = fp[i + 1, j]

        # getting the Fi value using the formula Fi = A^-1 * Fh
        Fic = inverse_A @ Fhc
        Fip = inverse_A @ Fhp

        # putting back the value of Fi back fo fc array
        fc[i] = Fic.T
        fp[i] = Fip.T

    # get the value of k
    k = int(np.floor(S / ds))

    # option price Call
    Vc = fc[0, k] + (fc[0, k + 1] - fc[0, k]) / ds * (S - K * ds)
    Vp = fp[0, k] + (fp[0, k + 1] - fp[0, k]) / ds * (S - K * ds)

    print('Call:', Vc, 'Put:', Vp)

    # Returning a tuple of call and put price
    return (Vc, Vp)

