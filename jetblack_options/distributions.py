"""Distributions"""

from math import exp, log, pi, sqrt

                                                                                         
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
            a *= y + 6.37396220353165
            a *= y + 33.912866078383
            a *= y + 112.079291497871
            a *= y + 221.213596169931
            a *= y + 220.206867912376
            b = 8.83883476483184E-02 * y + 1.75566716318264
            b *= y + 16.064177579207
            b *= y + 86.7807322029461
            b *= y + 296.564248779674
            b *= y + 637.333633378831
            b *= y + 793.826512519948
            b *= y + 440.413735824752
            c = e * a / b
        else:
            a = y + 0.65
            a = y + 4 / a
            a = y + 3 / a
            a = y + 2 / a
            a = y + 1 / a
            c = e / (a * 2.506628274631)
  
    return c - 1 if x > 0 else c


# Inverse cummulative normal distribution function
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
