Option Explicit     'Requires that all variables to be declared explicitly.
Option Compare Text 'Uppercase letters to be equivalent to lowercase letters.
Global Const Pi = 3.14159265358979


Option Base 1       'The "Option Base" statement allows to specify 0 or 1 as the
                    'default first index of arrays.

' Implementation By Espen Gaarder Haug
' Copyright Espen Gaarder Haug 2006





' What asset price that gives maximum DdeltaDvol
Public Function MaxDdeltaDvolAsset(UpperLowerFlag As String, x As Double, T As Double, b As Double, v As Double) As Double
    ' UpperLowerFlag"l" gives lower asset level that gives max DdeltaDvol
    ' UpperLowerFlag"l" gives upper asset level that gives max DdeltaDvol
    
    If UpperLowerFlag = "l" Then
        MaxDdeltaDvolAsset = x * Exp(-b * T - v * Sqr(T) * Sqr(4 + T * v ^ 2) / 2)
    ElseIf UpperLowerFlag = "u" Then
        MaxDdeltaDvolAsset = x * Exp(-b * T + v * Sqr(T) * Sqr(4 + T * v ^ 2) / 2)
    End If
    
End Function

' What strike price that gives maximum DdeltaDvol
Public Function MaxDdeltaDvolStrike(UpperLowerFlag As String, S As Double, T As Double, b As Double, v As Double) As Double
    
    ' UpperLowerFlag"l" gives lower strike level that gives max DdeltaDvol
    ' UpperLowerFlag"l" gives upper strike level that gives max DdeltaDvol

    If UpperLowerFlag = "l" Then
        MaxDdeltaDvolStrike = S * Exp(b * T - v * Sqr(T) * Sqr(4 + T * v ^ 2) / 2)
    ElseIf UpperLowerFlag = "u" Then
        MaxDdeltaDvolStrike = S * Exp(b * T + v * Sqr(T) * Sqr(4 + T * v ^ 2) / 2)
    End If
    
End Function

' What strike price that gives maximum gamma and vega
Public Function GMaxGammaVegaatX(S As Double, b As Double, T As Double, v As Double)

    GMaxGammaVegaatX = S * Exp((b + v * v / 2) * T)

End Function

' What asset price that gives maximum gamma
Public Function GMaxGammaatS(x As Double, b As Double, T As Double, v As Double)

    GMaxGammaatS = x * Exp((-b - 3 * v * v / 2) * T)

End Function

' What asset price that gives maximum vega
Public Function GMaxVegaatS(x As Double, b As Double, T As Double, v As Double)

    GMaxVegaatS = x * Exp((-b + v * v / 2) * T)
            
End Function

' Delta for the generalized Black and Scholes formula
Public Function GInTheMoneyProbability(CallPutFlag As String, S As Double, x As Double, T As Double, b As Double, v As Double) As Double

    Dim d2 As Double
    
    d2 = (Log(S / x) + (b - v ^ 2 / 2) * T) / (v * Sqr(T))
    
    If CallPutFlag = "c" Then
        GInTheMoneyProbability = CND(d2)
    ElseIf CallPutFlag = "p" Then
        GInTheMoneyProbability = CND(-d2)
    End If
    
End Function

' MirrorDeltaStrike, delta neutral straddle strike in the BSM formula
Public Function GDeltaMirrorStrike(S As Double, T As Double, b As Double, v As Double) As Double
    
    GDeltaMirrorStrike = S * Exp((b + v ^ 2 / 2) * T)
    
End Function

' MirrorProbabilityStrike, probability neutral straddle strike in the BSM formula
Public Function GProbabilityMirrorStrike(S As Double, T As Double, b As Double, v As Double) As Double

    GProbabilityMirrorStrike = S * Exp((b - v ^ 2 / 2) * T)
    
End Function

' MirrorDeltaStrike, general delta symmmetric strike in the BSM formula
Public Function GDeltaMirrorCallPutStrike(S As Double, x As Double, T As Double, b As Double, v As Double) As Double
    
    GDeltaMirrorCallPutStrike = S ^ 2 / x * Exp((2 * b + v ^ 2) * T)
    
End Function

' Elasticity for the generalized Black and Scholes formula
Public Function GElasticity(CallPutFlag As String, S As Double, x As Double, T As Double, r As Double, b As Double, v As Double) As Double

    GElasticity = GDelta(CallPutFlag, S, x, T, r, b, v) * S / GBlackScholes(CallPutFlag, S, x, T, r, b, v)
    
End Function

' Volatility estimate confidence interval
Function GConfidenceIntervalVolatility(Alfa As Double, n As Integer, VolatilityEstimate As Double, UpperLower As String)
    'UpperLower     ="L" gives the lower cofidence interval
    '               ="U" gives the upper cofidence interval
    'n: number of observations
    If UpperLower = "L" Then
        GConfidenceIntervalVolatility = VolatilityEstimate * Sqr((n - 1) / (Application.ChiInv(Alfa / 2, n - 1)))
    ElseIf UpperLower = "U" Then
        GConfidenceIntervalVolatility = VolatilityEstimate * Sqr((n - 1) / (Application.ChiInv(1 - Alfa / 2, n - 1)))
    End If

End Function


' Profitt/Loss STD for the generalized Black and Scholes formula
Public Function GProfitLossSTD(TypeFlag As String, CallPutFlag As String, S As Double, x As Double, T As Double, r As Double, b As Double, v As Double, NHedges As Integer) As Double
    
    If TypeFlag = "a" Then ' in dollars
        GProfitLossSTD = Sqr(Application.Pi() / 4) * GVega(S, x, T, r, b, v) * v / Sqr(NHedges)
    ElseIf TypeFlag = "p" Then ' in percent
        GProfitLossSTD = Sqr(Application.Pi() / 4) * GVega(S, x, T, r, b, v) * v / Sqr(NHedges) / GBlackScholes(CallPutFlag, S, x, T, r, b, v)
    End If

End Function

' GGammaPDtime for the generalized Black and Scholes formula
Public Function GDgammaPDtime(S As Double, x As Double, T As Double, r As Double, b As Double, v As Double) As Double
    
    Dim d1 As Double, d2 As Double
    
    d1 = (Log(S / x) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
    GDgammaPDtime = GGammaP(S, x, T, r, b, v) * (r - b + b * d1 / (v * Sqr(T)) + (1 - d1 * d2) / (2 * T))

End Function

Function GBlackScholesVarianceNGreeks(OutPutFlag As String, CallPutFlag As String, S As Double, x As Double, T As Double, r As Double, b As Double, v As Double, Optional dS)
            
    If IsMissing(dS) Then
        dS = 0.01
    End If
    
    
    If OutPutFlag = "p" Then ' Value
        GBlackScholesVarianceNGreeks = GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v)
    ElseIf OutPutFlag = "d" Then 'Delta
         GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S + dS, x, T, r, b, v) - GBlackScholesVariance(CallPutFlag, S - dS, x, T, r, b, v)) / (2 * dS)
    ElseIf OutPutFlag = "e" Then 'Elasticity
         GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S + dS, x, T, r, b, v) - GBlackScholesVariance(CallPutFlag, S - dS, x, T, r, b, v)) / (2 * dS) * S / GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v)
    ElseIf OutPutFlag = "g" Then 'Gamma
        GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S + dS, x, T, r, b, v) - 2 * GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v) + GBlackScholesVariance(CallPutFlag, S - dS, x, T, r, b, v)) / dS ^ 2
    ElseIf OutPutFlag = "gv" Then 'DGammaDvariance
        GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S + dS, x, T, r, b, v + 0.01) - 2 * GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v + 0.01) + GBlackScholesVariance(CallPutFlag, S - dS, x, T, r, b, v + 0.01) - GBlackScholesVariance(CallPutFlag, S + dS, x, T, r, b, v - 0.01) + 2 * GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v - 0.01) - GBlackScholesVariance(CallPutFlag, S - dS, x, T, r, b, v - 0.01)) / (2 * 0.01 * dS ^ 2) / 100
    ElseIf OutPutFlag = "gp" Then 'GammaP
        GBlackScholesVarianceNGreeks = S / 100 * (GBlackScholesVariance(CallPutFlag, S + dS, x, T, r, b, v) - 2 * GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v) + GBlackScholesVariance(CallPutFlag, S - dS, x, T, r, b, v)) / dS ^ 2
    ElseIf OutPutFlag = "dddv" Then 'DDeltaDvariance
        GBlackScholesVarianceNGreeks = 1 / (4 * dS * 0.01) * (GBlackScholesVariance(CallPutFlag, S + dS, x, T, r, b, v + 0.01) - GBlackScholesVariance(CallPutFlag, S + dS, x, T, r, b, v - 0.01) - GBlackScholesVariance(CallPutFlag, S - dS, x, T, r, b, v + 0.01) + GBlackScholesVariance(CallPutFlag, S - dS, x, T, r, b, v - 0.01)) / 100
    ElseIf OutPutFlag = "v" Then 'Variance Vega
        GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v + 0.01) - GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v - 0.01)) / 2
    ElseIf OutPutFlag = "vp" Then 'Variance VegaP
        GBlackScholesVarianceNGreeks = v / 0.1 * (GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v + 0.01) - GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v - 0.01)) / 2
    ElseIf OutPutFlag = "dvdv" Then 'Variance Dvegavariance
        GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v + 0.01) - 2 * GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v) + GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v - 0.01))
    ElseIf OutPutFlag = "t" Then 'Theta
        If T <= 1 / 365 Then
            GBlackScholesVarianceNGreeks = GBlackScholesVariance(CallPutFlag, S, x, 0.00001, r, b, v) - GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v)
        Else
            GBlackScholesVarianceNGreeks = GBlackScholesVariance(CallPutFlag, S, x, T - 1 / 365, r, b, v) - GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v)
        End If
    ElseIf OutPutFlag = "r" Then 'Rho
        GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S, x, T, r + 0.01, b + 0.01, v) - GBlackScholesVariance(CallPutFlag, S, x, T, r - 0.01, b - 0.01, v)) / (2)
    ElseIf OutPutFlag = "fr" Then 'Futures options rho
        GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S, x, T, r + 0.01, 0, v) - GBlackScholesVariance(CallPutFlag, S, x, T, r - 0.01, 0, v)) / (2)
    ElseIf OutPutFlag = "f" Then 'Rho2
        GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S, x, T, r, b - 0.01, v) - GBlackScholesVariance(CallPutFlag, S, x, T, r, b + 0.01, v)) / (2)
    ElseIf OutPutFlag = "b" Then 'Carry
        GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S, x, T, r, b + 0.01, v) - GBlackScholesVariance(CallPutFlag, S, x, T, r, b - 0.01, v)) / (2)
    ElseIf OutPutFlag = "s" Then 'Speed
        GBlackScholesVarianceNGreeks = 1 / dS ^ 3 * (GBlackScholesVariance(CallPutFlag, S + 2 * dS, x, T, r, b, v) - 3 * GBlackScholesVariance(CallPutFlag, S + dS, x, T, r, b, v) + 3 * GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v) - GBlackScholesVariance(CallPutFlag, S - dS, x, T, r, b, v))
    ElseIf OutPutFlag = "dx" Then 'Strike Delta
        GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S, x + dS, T, r, b, v) - GBlackScholesVariance(CallPutFlag, S, x - dS, T, r, b, v)) / (2 * dS)
    ElseIf OutPutFlag = "dxdx" Then 'Gamma
        GBlackScholesVarianceNGreeks = (GBlackScholesVariance(CallPutFlag, S, x + dS, T, r, b, v) - 2 * GBlackScholesVariance(CallPutFlag, S, x, T, r, b, v) + GBlackScholesVariance(CallPutFlag, S, x - dS, T, r, b, v)) / dS ^ 2
    End If

End Function

' The generalized Black and Scholes formula on variance form
Public Function GBlackScholesVariance(CallPutFlag As String, S As Double, x As Double, T As Double, r As Double, b As Double, v As Double) As Double

    Dim d1 As Double, d2 As Double
    d1 = (Log(S / x) + (b + v / 2) * T) / Sqr(v * T)
    d2 = d1 - Sqr(v * T)

    If CallPutFlag = "c" Then
        GBlackScholesVariance = S * Exp((b - r) * T) * CND(d1) - x * Exp(-r * T) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        GBlackScholesVariance = x * Exp(-r * T) * CND(-d2) - S * Exp((b - r) * T) * CND(-d1)
    End If
    
End Function
