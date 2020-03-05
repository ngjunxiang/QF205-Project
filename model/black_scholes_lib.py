
def check_params(S, K, r, q, T, sigma, M, N):
    if S >= 2 * K or N == 0 or M == 0:
        return False

    return True