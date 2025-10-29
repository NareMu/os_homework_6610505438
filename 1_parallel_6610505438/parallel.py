"""
Parallel of Quadratic Sieve (QS) 
write from idea of this page https://www.bytopia.dk/qs/
"""
import math
import time
import sys
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional
from mpi4py import MPI
import sympy.ntheory as nt
from sympy.functions.combinatorial.numbers import legendre_symbol

# --------------------------------- Tags  ----------------------------------------
AFXTAG = 1
AXBTAG = 2
RESULTACKTAG = 4
DTAG = 8
NTAG = 12
REQDTAG = 16
REQNTAG = 20
RESULTTAG = 24
DIETAG = 28

EXTRA = 10
SIEVE_LENGTH = 500_000

# -------------------------- Math -------------------
def tonelli_shanks(n: int, p: int) -> int:
    """perform Tonelli-Shanks algorithm from given number(n) and prime(p)"""
    assert p > 2 ,f"tonelli_shanks expected odd prime but {p} given"
    n %= p
    if n == 0:
        return 0
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    # factor p-1 = q*2^s
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    # find a quadratic non-residue z
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)
    while t != 1:
        # find least i (0<i<m) with t^(2^i) == 1
        i, t2i = 1, (t * t) % p
        while i < m and t2i != 1:
            t2i = (t2i * t2i) % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r

# -------------------------- Factor base & polynomial data -----------------------

@dataclass
class FBElement:
    p: int
    logp: float
    sqrtn1: int
    sqrtn2: int
    y: int = 0  
    z: int = 0  

@dataclass
class QSParams:
    N: int
    factor_base_size: int
    T: float
    M: int

def init_params(number_str: str) -> QSParams:
    """ initialize a needed params from given number in string and return in QSParams"""
    x = len(number_str)
    if x < 9:
        raise ValueError("Number too small for QS")
    N = int(number_str)
    if x < 25:
        fbs = 100
    else:
        fbs = int(2.93 * (x * x) - 164.4 * x + 2455)
    M = int(386 * (x * x) - 23209.3 * x + 352768)
    T = 0.0268849 * x + 0.783929
    return QSParams(N=N, factor_base_size=fbs, T=T, M=M)

def compute_factor_base(params: QSParams) -> Tuple[List[FBElement], int]:
    """form a valid factor base of given number"""
    N = params.N
    factor_base_size = params.factor_base_size
    fb: List[FBElement] = [FBElement(p=1, logp=0, sqrtn1=0, sqrtn2=0)]  
    p = 3
    while len(fb) < factor_base_size:
        if legendre_symbol(N, p) == 1:
            nmodp = N % p
            r = tonelli_shanks(nmodp, p)
            fb.append(FBElement(p=p, logp=int(math.log(p, 2)), sqrtn1=r, sqrtn2=(p - r) ))
        p = nt.nextprime(p)
    # threshold log(M*sqrt(N/2) / p_max^T) in bits
    tmp = params.M * int(math.isqrt(math.ceil(N/2)))
    cp = params.T * (math.log(fb[-1].p, 2))
    pmaxT = 2 ** int(cp)
    threshold_bits = (math.ceil(tmp/pmaxT)).bit_length()
    return fb, threshold_bits

# -------------------------- Polynomial selection --------------------------------

def get_next_D(N: int, prev_D: int) -> int:
    """find next valid D that use for form a polynomial from a given N and previous D """
    while True:
        D = nt.nextprime(prev_D)
        if D % 4 == 3 and legendre_symbol(N, D) == 1:
            return D

def generate_polynomial(N: int, D: int) -> Tuple[int, int, int]:
    """find A,B,C that use for form a polynomial from a given N and D"""
    A = D * D
    h1 = pow(N, (D + 1) // 4, D)
    tmp = (N - (h1 * h1)) // D
    inv_2h1 = pow((2 * h1) , -1, D)
    h2 = (inv_2h1 * (tmp % D)) % D
    B = (h1 + h2 * D) % A
    C =  (B * B - N) // A
    return A, B, C

def initialize_polynomial(A: int, B: int, fb: List[FBElement]) -> None:
    """find valid root of each factor base"""
    for i in range(1, len(fb)):
        p = fb[i].p
        if (2 * A) % p == 0:
            fb[i].y = fb[i].z = 0
            continue
        inv_down = pow(2 * A , -1, p)
        # find two root from quadratic equation
        up1 = (-2 * B + 2 * fb[i].sqrtn1) % p
        fb[i].y = (up1 * inv_down) % p
        up2 = (-2 * B + 2 * fb[i].sqrtn2) % p
        fb[i].z = (up2 * inv_down) % p

# -------------------------- Sieve & trial division --------------------------------

def poly_val(x: int, A: int, B: int, C: int) -> int:
    """calculate a value of polynomial from a given x, A, B, C"""
    return A * x * x + 2 * B * x + C

def sieve_interval(A: int, B: int, C: int, fb: List[FBElement], threshold_bits: float,
                   start_offset: int, length: int) -> List[Tuple[int, int, int]]:
    """Find valid candidate for f(x) then return in a list of tuple that contain(vec,axb,afx) which 
    vec - a vector of power of prime factor in GF(2)
    axb - (ax+b)^2 
    afx- Af(x)
    """
    logfxs = np.zeros(length, dtype=float)
    # precompute in log to cutoff some impossible fx for reduce large number dividing
    for i in range(1, len(fb)):
        p = fb[i].p
        logp = fb[i].logp
        y, z = fb[i].y, fb[i].z
        if p <= 0:
            continue
        q = (y - start_offset) % p
        if q < length:
            logfxs[q: length: p] += logp
        q = (z - start_offset) % p
        if q < length:
            logfxs[q: length: p] += logp

    hits: List[Tuple[int, int, int]] = []
    # mark a number of each factor base contain in fx and also mark right most bit to be a sign bit in GF(2)
    for i in range(length):
        if logfxs[i] >= threshold_bits:
            x = start_offset + i
            fx = poly_val(x, A, B, C)
            vec = 0
            val = fx
            if val < 0:
                vec = 1
                val = -val
            for j in range(1, len(fb)):
                p = fb[j].p
                while val % p == 0:
                    vec ^= (1 << j)
                    val //= p
            if val == 1:
                axb = A * x + B
                a_fx = A * fx
                hits.append((vec, axb, a_fx))
    return hits


# -------------------------- Gaussian elimination over GF(2) ----------------------

def gauss_eliminate(ncols: int, rows: List[int]) -> List[int]:
    """perform gaussian elimination in GF(2) from a given metrix via bitwise operation and return a result metrix in form of list of int (each int refer to a list of bit)"""
    augmented = []
    nrow = len(rows)
    # form a augmented matrix [A|I] 
    for i, r in enumerate(rows):
        comb = 1 << i
        augmented.append(r << nrow | comb)

    col = 0
    # gauss_elimination
    for r in range(nrow):
        if col >= ncols:
            break
        pivot = None
        for i in range(r, len(augmented)):
            if (augmented[i] >> (col + nrow)) & 1:
                pivot = i
                break
        if pivot is None:
            col += 1
            continue
        if pivot != r:
            augmented[r], augmented[pivot] = augmented[pivot], augmented[r]

        for i in range(r + 1, len(augmented)):
            if (augmented[i] >> (col + nrow)) & 1:
                augmented[i] ^= augmented[r]
        col += 1
    return augmented

def back_tracking(axb: np.ndarray, afx: np.ndarray, augmented: np.ndarray, ncol: int, nrow: int, N:int) -> Optional[Tuple[int, int]]:
    """use a result from gaussian elimination to find a valid answers"""
    # masking for get only A part from [A|I]
    A_mask = ((1 << (ncol-nrow)) - 1) << nrow
    for row in augmented:
        # if A is not all zero then this row can't be the answer
        if (row & A_mask) != 0:
            continue 
        x = 1
        y = 1
        # choose a row j when in I column j is 1
        for j in range(ncol):
            if (row >> j) & 1:
                x *= axb[j]
                y *= afx[j]

        y = int(math.isqrt(y))
        # test it is a valid answer
        g = math.gcd(y - x, N)
        if 1 < g < N:
            p = g
            q = N // g
            return (p,q)
    return None


# -------------------------- Master/Slave protocol --------------------------------
def master_main_serial(argv: List[str], comm: MPI.Comm) -> None:
    """perform a qs algorithm with non-parallel method"""
    rank = comm.Get_rank()
    assert rank == 0, f"rank: {rank} can't be master"

    if len(argv) != 2:
        print("Usage: parallel.py <number to factor>")
        return

    params = init_params(argv[1])
    N = params.N
    M = params.M
    factor_base_size = params.factor_base_size
    if nt.isprime(N):
        print(f"{N} is prime")
        return

    fb, threshold_bits = compute_factor_base(params)

    print("Sieving without slave nodes ...")
    print(f"Factor base size: {factor_base_size}, threshold(bits): {threshold_bits:.2f}",flush = True)

    needed = factor_base_size + EXTRA

    exponent_matrix = np.zeros(needed, dtype=object)  
    axb = np.zeros(needed, dtype=object)              
    afx = np.zeros(needed, dtype=object)              

    next_index = 0
    D_prev = max(int(math.isqrt(int(math.isqrt(2 * N)/M)))-1,3)

    t0 = time.time()

    while next_index < needed:
        D = get_next_D(N, D_prev)
        D_prev = D

        A, B, C = generate_polynomial(N, D)
        initialize_polynomial(A, B, fb)

        start = -params.M
        remaining = 2 * params.M

        while remaining > 0 and next_index < needed:
            length = min(SIEVE_LENGTH, remaining)
            hits = sieve_interval(A, B, C, fb, threshold_bits, start, length)

            for vec, axb_val, a_fx_val in hits:
                if next_index >= needed:
                    break
                exponent_matrix[next_index] = vec
                axb[next_index] = axb_val
                afx[next_index] = a_fx_val
                next_index += 1

            start += length
            remaining -= length

    print()
    t_sieve = time.time() - t0
    print(f"Done sieving ({t_sieve} s)\n")

    print("Performing Gauss elimination ...")
    t1 = time.time()

    # perform a gauss elimination to find some combination that prod(Af(xi)) can write in form prod of factor_base with even power
    augmented = gauss_eliminate(factor_base_size,list(exponent_matrix[:needed]))
    t2 = time.time()
    print(f"Done ({t2 - t1} s)\n")

    nrow = needed
    ncol = factor_base_size + nrow
    result = back_tracking(axb, afx, augmented, ncol, nrow, N)

    if result is None:
        print("Unable to find a non-trivial factorization from dependencies. Try more relations or adjust params.")
    else:
        p,q = result
        print("The factorization is:\n")
        print(f"{N} = {p} * {q}\n")
    print(f"Total running time: {time.time() - t0} s")

def master_main(argv: List[str], comm: MPI.Comm) -> None:
    """perform a qs algorithm with parallel method, this function will perform as a master(process0)"""
    size = comm.Get_size()
    rank = comm.Get_rank()
    assert rank == 0 ,f"rank: {rank} can't be master"

    if len(argv) != 2:
        print("Usage: parallel.py <number to factor>")
        return

    params = init_params(argv[1])
    N = params.N
    M = params.M
    factor_base_size = params.factor_base_size
    if nt.isprime(N):
        print(f"{N} is prime")
        return
    _, threshold_bits = compute_factor_base(params)

    print(f"Sieving with {size-1} slave node ...")
    print(f"  Factor base size: {factor_base_size}, threshold(bits): {threshold_bits:.2f}", flush=True)

    relations = 0
    next_index = 0
    needed = factor_base_size + EXTRA
    D_prev = max(int(math.isqrt(int(math.isqrt(2 * N)/M)))-1,3)
    
    exponent_matrix = np.zeros(needed, dtype=object)
    axb = np.zeros(needed, dtype=object)
    afx = np.zeros(needed, dtype=object)

    terminated = 0
    t0 = time.time()

    while relations < needed or terminated < size - 1:
        status = MPI.Status()
        msg = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()
        src = status.Get_source()
        # slave request N
        if tag == REQNTAG:
            if next_index >= needed:
                comm.send(None, dest=src, tag=DIETAG)
                terminated += 1
            else:
                comm.send(N, dest=src, tag=NTAG)
        # slave request D
        elif tag == REQDTAG:
            if next_index >= needed:
                comm.send(None, dest=src, tag=DIETAG)
                terminated += 1
            else:
                D = get_next_D(N,D_prev)
                D_prev = D
                comm.send(D, dest=src, tag=DTAG)
        # slave request index
        elif tag == RESULTTAG:
            if next_index >= needed:
                comm.send(-1, dest=src, tag=DIETAG)
                terminated += 1
            else:
                vec = msg
                exponent_matrix[next_index] = vec
                comm.send(next_index, dest=src, tag=RESULTACKTAG)
                next_index += 1
        # slave send an answer
        else:
            client_index = tag >> 2
            which = tag & 3
            if which == AXBTAG:
                axb[client_index] = msg
            elif which == AFXTAG:
                afx[client_index] = msg
                relations += 1

    print()
    t_sieve = time.time() - t0
    print(f"Done sieving ({t_sieve} s)\n")

    print("Performing Gauss elimination ...")
    t1 = time.time()

    # perform a gauss elimination to find some combination that prod(Af(x)) is can write in form prod of factor_base with even power
    augmented = gauss_eliminate(factor_base_size, list(exponent_matrix[:needed]))
    t2 = time.time()
    print(f"Done ({t2 - t1} s)\n")

    nrow = needed
    ncol = factor_base_size + nrow
    result = back_tracking(axb, afx, augmented, ncol, nrow, N)

    if result is None:
        print("Unable to find a non-trivial factorization from dependencies. Try more relations or adjust params.")
    else:
        p,q = result
        print("The factorization is:\n")
        print(f"{N} = {p} * {q}\n")

    print(f"Total running time: {time.time() - t0} s")

def slave_main(comm: MPI.Comm) -> None:
    """perform a qs algorithm with parallel method, this function will perform as a slave(all process exept process 0)"""
    rank = comm.Get_rank()

    print(f"rank: {rank} start process",flush = True)
    comm.send(None, dest=0, tag=REQNTAG)
    status = MPI.Status()
    N = comm.recv(source=0, tag=NTAG)

    params = init_params(str(N))
    fb, threshold_bits = compute_factor_base(params)

    comm.send(None, dest=0, tag=REQDTAG)
    D = comm.recv(source=0, tag=DTAG)

    while True:
        A, B, C = generate_polynomial(params.N, D)
        initialize_polynomial(A, B, fb)

        start = -params.M
        remaining = 2 * params.M
        while remaining > 0:
            length = min(SIEVE_LENGTH, remaining)
            hits = sieve_interval(A, B, C, fb, threshold_bits, start, length)
            for vec, axb_val, a_fx_val in hits:
                # send request the next index
                comm.send(vec, dest=0, tag=RESULTTAG)
                status = MPI.Status()
                # receive index
                idx = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
                if status.Get_tag() != RESULTACKTAG:
                    return
                # send real data with index
                comm.send(axb_val, dest=0, tag=((idx << 2) + AXBTAG))
                comm.send(a_fx_val, dest=0, tag=((idx << 2) + AFXTAG))
            start += length
            remaining -= length

        comm.send(None, dest=0, tag=REQDTAG)
        status = MPI.Status()
        msg = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        if status.Get_tag() != DTAG:
            return
        D = msg

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # assign work to each process following their rank and total size of process
    if rank == 0 and size > 1:
        master_main(sys.argv, comm)
    elif rank == 0:
        master_main_serial(sys.argv,comm)
    else:
        slave_main(comm)
    MPI.Finalize()


if __name__ == "__main__":
    main()

