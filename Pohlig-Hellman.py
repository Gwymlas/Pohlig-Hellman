import gmpy2


# def gcd(a, b):
#     while b != 0:
#         a, b = b, a % b
#     return a


def multiplicative_order(g, p):
    assert gmpy2.gcd(g, p) == 1

    result = gmpy2.mpz(1)
    order = gmpy2.mpz(1)
    while order < p:
        # result = (result * g) % p
        result = gmpy2.mod(gmpy2.mul(result, g), p)

        if result == 1:
            return order

        order += 1

    return -1


def multiplicative_inv(g, p):
    d, inv, m = gmpy2.gcdext(g, p)
    assert d == 1
    return gmpy2.mod(inv, p)


def chinese_remainder_theorem(y, m):
    assert len(y) == len(m)

    x = gmpy2.mpz(0)
    M = gmpy2.mpz(1)

    for m_i in m:
        M = gmpy2.mul(M, m_i)

    for i in range(len(y)):
        b = gmpy2.mpz(gmpy2.div(M, m[i]))
        c = multiplicative_inv(b, m[i])
        x += gmpy2.mod(gmpy2.mul(gmpy2.mul(y[i], b), c), M)

    return x, M


def brute_force(g, h, p):
    value = gmpy2.mpz(1)
    result = gmpy2.mpz(0)

    if gmpy2.is_prime(p):
        order = p - 1
    else:
        order = multiplicative_order(g, p)

    assert order != -1

    while result < order:
        value = gmpy2.mod(gmpy2.mul(value, g), p)
        result += 1

        if value == h:
            return result, order

    return None


def bsgs(g, h, p):
    if gmpy2.is_prime(p):
        order = p - 1
    else:
        order = multiplicative_order(g, p)
        assert order != -1

    n = 1 + gmpy2.isqrt(order)

    # работает в 5.5 раза дольше
    # baby_step = {pow(g, i, mod=p): i for i in range(n)}

    baby_step = dict()
    elem_baby_step = gmpy2.mpz(1)

    for i in range(n):
        baby_step[elem_baby_step] = i
        elem_baby_step = gmpy2.mod(gmpy2.mul(elem_baby_step, g), p)

    multiplier = gmpy2.powmod(g, -n, p)
    elem = h

    for j in range(n):
        if elem in baby_step:
            return (j * n + baby_step[elem]) % p, order

        elem = gmpy2.mod(gmpy2.mul(elem, multiplier), p)

    return None, order


def pohlig_hellman_prime(g, h, p, q, n):
    order = multiplicative_order(g, p)
    assert order == gmpy2.powmod(q, n, p)

    g_i = [gmpy2.mpz(0) for _ in range(n)]
    h_i = [gmpy2.mpz(0) for _ in range(n)]

    g_i[n-1] = g
    h_i[n-1] = h

    for i in range(1, n):
        g_i[n - 1 - i] = gmpy2.powmod(g_i[n-i], q, p)
        h_i[n - 1 - i] = gmpy2.powmod(h_i[n - i], q, p)

    x_i = [gmpy2.mpz(0) for _ in range(n)]
    x_i[0], _ = brute_force(g_i[0], h_i[0], p)

    for i in range(1, n):
        # gamma = math.prod([pow(g_i[i-j], x_i[j], p) for j in range(i)]) % p
        gamma = gmpy2.mpz(1)
        for j in range(i):
            gamma = gmpy2.mod(gmpy2.mul(gamma, gmpy2.powmod(g_i[i-j], x_i[j], p)), p)

        beta = gmpy2.mod(gmpy2.mul(h_i[i], multiplicative_inv(gamma, p)), p)
        x_i[i], _ = brute_force(g_i[0], beta, p)

    x = 0

    for i, elem in enumerate(x_i):
        x += gmpy2.mul(elem, gmpy2.powmod(q, i, p))

    return gmpy2.mod(x, order)


def pohlig_hellman(g, h, p, q, n):
    assert len(q) == len(n)
    order = multiplicative_order(g, p)

    g_i = [gmpy2.mpz(0) for _ in range(len(q))]
    h_i = [gmpy2.mpz(0) for _ in range(len(q))]
    m_i = [gmpy2.mpz(0) for _ in range(len(q))]

    for i in range(len(q)):
        m_i[i] = gmpy2.powmod(q[i], n[i], p)
        power = gmpy2.mpz(gmpy2.div(order, m_i[i]))
        g_i[i] = gmpy2.powmod(g, power, p)
        h_i[i] = gmpy2.powmod(h, power, p)

    y_i = [pohlig_hellman_prime(g_i[i], h_i[i], p, q[i], n[i]) for i in range(len(q))]

    x, mod = chinese_remainder_theorem(y_i, m_i)
    return x, mod


def main():
    g = gmpy2.mpz("2")
    h = gmpy2.mpz("52")
    p = gmpy2.mpz("61")

    # print(multiplicative_order(g, p))

    # print(multiplicative_order(g, p))

    # print(brute_force(g, h, p))

    # print(bsgs(g, h, p))

    print(pohlig_hellman_prime(2, 193, 257, 2, 4))

    print(pohlig_hellman(g, h, p,
                         [gmpy2.mpz(2), gmpy2.mpz(3), gmpy2.mpz(5)],
                         [gmpy2.mpz(2), gmpy2.mpz(1), gmpy2.mpz(1)]))


if __name__ == '__main__':
    main()
