from sage.all import *
import time

def coppersmith_howgrave_univariate(pol, modulus, beta, mm, tt, XX, debug=False):
    """
    Implements Coppersmith's method with Howgrave-Graham's optimization.
    
    :param pol: Polynomial f(x) we want to find small roots of.
    :param modulus: RSA modulus (N)
    :param beta: Parameter β where 0 < β ≤ 1
    :param mm: Lattice row parameter (optimized based on β)
    :param tt: Number of extra polynomials used
    :param XX: Upper bound for the small root
    :param debug: If True, display intermediate computations
    :return: List of small roots found
    """
    dd = pol.degree()
    nn = dd * mm + tt  # Total rows in lattice

    if not pol.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    # Construct polynomials
    polZ = pol.change_ring(ZZ)  # Convert to integer coefficients
    x = polZ.parent().gen()

    gg = []
    for i in range(mm):
        for j in range(dd):
            gg.append((x * XX)**j * modulus**(mm - i) * polZ(x * XX)**i)
    for i in range(tt):
        gg.append((x * XX)**i * polZ(x * XX)**mm)

    # Construct the lattice
    BB = Matrix(ZZ, nn)
    for i in range(nn):
        for j in range(i + 1):
            BB[i, j] = gg[i][j]

    # Apply LLL reduction (without printing the matrix)
    BB = BB.LLL()

    # Construct new polynomial from smallest vector
    new_pol = sum(x**i * BB[0, i] / XX**i for i in range(nn))

    # Find integer roots
    potential_roots = new_pol.roots()
    roots = [ZZ(root[0]) for root in potential_roots if root[0].is_integer()]
    
    return roots


def attack_partial_prime(N, partial_p, unknown_bits, beta=0.5, debug=False):
    """
    Coppersmith's attack on RSA when a part of a prime factor is known.
    
    :param N: RSA modulus (N = p * q)
    :param partial_p: Known most significant bits of p
    :param unknown_bits: Number of bits missing from p
    :param beta: Parameter β where 0 < β ≤ 1
    :param debug: If True, display intermediate computations
    :return: Recovered prime factor p
    """
    F.<x> = PolynomialRing(Zmod(N))
    pol = x + partial_p  # f(x) = x + p0 (where x is the unknown part)
    dd = pol.degree()

    # Compute optimized parameters
    epsilon = beta / 7
    mm = ceil(beta**2 / (dd * epsilon))
    tt = floor(dd * mm * ((1/beta) - 1))
    XX = ceil(N**((beta**2/dd) - epsilon))

    if debug:
        print(f"Optimized parameters: m={mm}, t={tt}, X={XX}")

    # Run Coppersmith's method
    roots = coppersmith_howgrave_univariate(pol, N, beta, mm, tt, XX, debug)

    for root in roots:
        p_candidate = partial_p + root
        if N % p_candidate == 0:
            return p_candidate
    return None


# ======= TEST CASE 1: FACTORING WITH PARTIAL P KNOWN =======
print("\n[TEST 1] Recovering p with missing bits\n")

bits = 512
p = random_prime(2**bits)
q = random_prime(2**bits)
N = p * q

leak_bits = 50  # Assume we know the most significant bits of p
partial_p = p - (p % 2**leak_bits)  # Remove last `leak_bits` bits

print(f"Original p:  {p}")
print(f"Known part:  {partial_p}")

# Run the attack
start_time = time.time()
recovered_p = attack_partial_prime(N, partial_p, leak_bits, beta=0.5, debug=False)

if recovered_p:
    print(f"Recovered p: {recovered_p}")
    assert recovered_p == p, "Recovered p is incorrect!"
    print("✅ Attack successful!")
else:
    print("❌ Attack failed.")

print(f"Time taken: {time.time() - start_time:.2f} seconds")


# ======= TEST CASE 2: BREAKING LOW EXPONENT RSA =======
print("\n[TEST 2] Recovering small e RSA message\n")

N = random_prime(2**512) * random_prime(2**512)
e = 3  # Small exponent
M = ZZ.random_element(2**100)  # Random message
C = pow(M, e, N)

F.<x> = PolynomialRing(Zmod(N))
pol = x**e - C

# Run the attack
roots = coppersmith_howgrave_univariate(pol, N, beta=1, mm=8, tt=5, XX=2**100, debug=False)

if roots and roots[0] == M:
    print(f"Recovered message: {roots[0]}")
    print("✅ Attack successful!")
else:
    print("❌ Attack failed.")


# ======= TEST CASE 3: RECOVERING PARTIAL MESSAGE FROM STRUCTURED PADDING =======
print("\n[TEST 3] Recovering hidden message from padding\n")

N = random_prime(2**512) * random_prime(2**512)
e = 3  # Small exponent

# Message has known prefix "CRYPTO" and hidden suffix
prefix = int.from_bytes(b"CRYPTO", "big") << 100  # Known 48-bit prefix
hidden = ZZ.random_element(2**100)  # 100-bit hidden part
M = prefix + hidden  # Construct full message
C = pow(M, e, N)  # Encrypt

F.<x> = PolynomialRing(Zmod(N))
pol = (prefix + x)**e - C

# Run the attack
roots = coppersmith_howgrave_univariate(pol, N, beta=1, mm=8, tt=5, XX=2**100, debug=False)

if roots and prefix + roots[0] == M:
    print(f"Recovered hidden part: {roots[0]}")
    print("✅ Attack successful!")
else:
    print("❌ Attack failed.")
