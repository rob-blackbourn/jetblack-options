"""Distributions"""

from math import asin, exp, log, nan, pi, sin, sqrt
from statistics import NormalDist

NORMAL_DIST = NormalDist()
PDF = NORMAL_DIST.pdf
CDF = NORMAL_DIST.cdf
INV_CDF = NORMAL_DIST.inv_cdf
                                                                                         
# The normal distribution function
def ND(x: float) -> float:
    return 1 / sqrt(2 * pi) * exp(-x ** 2 / 2)


# Cumulative double precision algorithm based on Hart 1968
# Based on implementation by Graeme West
def CND(x: float)-> float:
    
    y = abs(x)
    if y > 37:
        return 0
    else:
        e = exp(-y ** 2 / 2)
        if y < 7.07106781186547:
            a = 3.52624965998911E-02 * y + 0.700383064443688
            a = a * y + 6.37396220353165
            a = a * y + 33.912866078383
            a = a * y + 112.079291497871
            a = a * y + 221.213596169931
            a = a * y + 220.206867912376
            b = 8.83883476483184E-02 * y + 1.75566716318264
            b = b * y + 16.064177579207
            b = b * y + 86.7807322029461
            b = b * y + 296.564248779674
            b = b * y + 637.333633378831
            b = b * y + 793.826512519948
            b = b * y + 440.413735824752
            c = e * a / b
        else:
            a = y + 0.65
            a = y + 4 / a
            a = y + 3 / a
            a = y + 2 / a
            a = y + 1 / a
            c = e / (a * 2.506628274631)
  
    return 1 - c if x > 0 else c


# Inverse cumulative normal distribution function
def CNDEV(U: float) -> float:
    
    A = (2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637)
    b = (-8.4735109309, 23.08336743743, -21.06224101826, 3.13082909833)
    c = (0.337475482272615, 0.976169019091719, 0.160797971491821, 2.76438810333863E-02, 3.8405729373609E-03, 3.951896511919E-04, 3.21767881767818E-05, 2.888167364E-07, 3.960315187E-07)

    x = U - 0.5
    if abs(x) < 0.92:
        r = x * x
        r = x * (((A[3] * r + A[2]) * r + A[1]) * r + A[0]) / ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1)
        return r

    r = U
    if x >= 0:
        r = 1 - U
    r = log(-log(r))
    r = c[0] + r * (c[1] + r * (c[2] + r * (c[3] + r + (c[4] + r * (c[5] + r * (c[6] + r * (c[7] + r * c[8])))))))
    if x < 0:
        r = -r
    return r

# Approximation of inverse chi-square (Weiss and Greenhall, 1996)
# Restrictions: df >= 1 and 0.005 <= p <= 0.995
# df can in theory be non-integral, but we define it as int here. 
# Returns, given the desired p-value and a degree of freedom df, 
# the corresponding chi-square value.
# Max. error of approximation is 3%
# See: https://apps.dtic.mil/sti/pdfs/ADA515532.pdf
def CHIINV(p: float, df: int) -> float:
    # double a, A, G, g, y, u, c1, c2, c3, c4, c5, x;
    # double p1, a0, a1, b1, b2, t, X, s, b ;
    # int k, n, i;
    p = 1.0 - p
    if p <= 0.5 and df <= 10:
        c1 = -0.5748646
        c2 =  0.9512363
        c3 = -0.6998588
        c4 =  0.4245549
        c5 = -0.1010678

        a = df / 2.0
        n = int(a)
        y = a - n
        G = 1 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * c5))))
        for k in range(1, n+1):
            G = G * (y + k)
        A = p * G
        u = 0.0
        for _ in range(8):
            g = 1 + (u/(a+1)) * ( 1 + (u/(a+2)) * (1 + (u/(a+3)) ))
            u = pow((A * exp(u)/g), 1/a)
        x = 2 * u
    else:
        a0 = 2.30753
        a1 = 0.27601
        b1 = 0.99229
        b2 = 0.04481
        
        p1 = min(p, 1.0 - p)
        t = sqrt(-2 * log(p1))
        X = t - (a0 + a1 * t) / (1 + b1 * t + b2 * t * t)
        s = -1 if p - 0.5 < 0 else 1 if p - 0.5 > 0 else 0
        #s = ( (p - 0.5) > 0 ) - ( (p-0.5) < 0 ) # sign(p-0.5)
        b = 2.0 / (9.0 * df)
        x = df * pow((1 - b + s * X * sqrt(b)), 3)

    return x

# The cumulative bivariate normal distribution function
def CBND(x: float, y: float, rho: float) -> float:
    #     A function for computing bivariate normal probabilities.
    #
    #       Alan Genz
    #       Department of Mathematics
    #       Washington State University
    #       Pullman, WA 99164-3113
    #       Email : alangenz@wsu.edu
    #
    #    This function is based on the method described by
    #        Drezner, Z and G.O. Wesolowsky, (1990),
    #        On the computation of the bivariate normal integral,
    #        Journal of Statist. Comput. Simul. 35, pp. 101-107,
    #    with major modifications for double precision, and for |R| close to 1.
    #   This code was originally translated into VBA by Graeme West

    if abs(rho) < 0.3:
        W = [
            0.17132449237917,
            0.360761573048138,
            0.46791393457269
        ]
        XX = [
            -0.932469514203152,
            -0.661209386466265,
            -0.238619186083197
        ]
    elif abs(rho) < 0.75:
        W = [
            4.71753363865118E-02,
            0.106939325995318,
            0.160078328543346,
            0.203167426723066,
            0.233492536538355,
            0.249147045813403
        ]
        XX = [
            -0.981560634246719,
            -0.904117256370475,
            -0.769902674194305,
            -0.587317954286617,
            -0.36783149899818,
            -0.125233408511469
        ]
    else:
        W = [
            1.76140071391521E-02,
            4.06014298003869E-02,
            6.26720483341091E-02,
            8.32767415767048E-02,
            0.10193011981724,
            0.118194531961518,
            0.131688638449177,
            0.142096109318382,
            0.149172986472604,
            0.152753387130726,
        ]
        XX = [
            -0.993128599185095,
            -0.963971927277914,
            -0.912234428251326,
            -0.839116971822219,
            -0.746331906460151,
            -0.636053680726515,
            -0.510867001950827,
            -0.37370608871542,
            -0.227785851141645,
            -7.65265211334973E-02
        ]
          
    h = -x
    k = -y
    hk = h * k
    BVN = 0.0
          
    if abs(rho) < 0.925:
        if abs(rho) > 0:
            hs = (h * h + k * k) / 2
            asr = asin(rho)
            for i in range(len(W)):
                for ISs in (-1, 1):
                    sn = sin(asr * (ISs * XX[i] + 1) / 2)
                    BVN = BVN + W[i] * exp((sn * hk - hs) / (1 - sn * sn))
            BVN = BVN * asr / (4 * pi)
        BVN = BVN + CND(-h) * CND(-k)
    else:
        if rho < 0:
            k = -k
            hk = -hk
        if abs(rho) < 1:
            Ass = (1 - rho) * (1 + rho)
            A = sqrt(Ass)
            bs = (h - k) ** 2
            c = (4 - hk) / 8
            d = (12 - hk) / 16
            asr = -(bs / Ass + hk) / 2
            if asr > -100:
                BVN = A * exp(asr) * (1 - c * (bs - Ass) * (1 - d * bs / 5) / 3 + c * d * Ass * Ass / 5)
            if -hk < 100:
                b = sqrt(bs)
                BVN = BVN - exp(-hk / 2) * sqrt(2 * pi) * CND(-b / A) * b * (1 - c * bs * (1 - d * bs / 5) / 3)
            A = A / 2
            for i in range(len(W)):
                for ISs in (-1, 1):
                    xs = (A * (ISs * XX[i] + 1)) ** 2
                    rs = sqrt(1 - xs)
                    asr = -(bs / xs + hk) / 2
                    if asr > -100:
                        BVN = BVN + A * W[i] * exp(asr) * (
                            exp(-hk * (1 - rs) / (2 * (1 + rs))) / rs -
                            (1 + c * xs * (1 + d * xs))
                        )
            BVN = -BVN / (2 * pi)
        if rho > 0:
            BVN = BVN + CND(-max(h, k))
        else:
            BVN = -BVN
            if k > h:
                BVN = BVN + CND(k) - CND(h)
    return BVN

