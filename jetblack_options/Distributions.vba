Global Const Pi = 3.14159265358979

Option Explicit     'Requires that all variables to be declared explicitly.
Option Compare Text 'Uppercase letters to be equivalent to lowercase letters.

Option Base 0       'The "Option Base" statment alowws to specify 0 or 1 as the
                             'default first index of arrays.

' Programmer Espen Gaarder Haug
' This code comes with The Complete Guide to Option Pricing Formulas, McGraw-Hill
' Copyright Espen G. Haug


' Convert rate to countinuous compounding rate
Public Function ConvertingToCCRate(r As Double, Compoundings As Double) As Double
    If Compoundings = 0 Then
        ConvertingToCCRate = r
    Else
        ConvertingToCCRate = Compoundings * Log(1 + r / Compoundings)
    End If
End Function
                    
' The cumulative bivariate normal distribution function
' Drezner-Wesolowsky 1990 simple algorithm
Public Function CBND2(A As Double, b As Double, rho As Double) As Double

    Dim g As Double, P As Double, x, y, sum As Double
    Dim i As Integer
    
    x = Array(0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042)
    y = Array(0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992)
    
    sum = 0
    For i = 0 To 4
        P = y(i) * rho
        g = 1 - P ^ 2
        sum = sum + x(i) * Exp((2 * A * b * P - A ^ 2 - b ^ 2) / g / 2) / Sqr(g)
    Next
    CBND2 = rho * sum + CND(A) * CND(b)
End Function

Public Function Max(x, y)
    Max = Application.Max(x, y)
End Function

Public Function Min(x, y)
    Min = Application.Min(x, y)
End Function




' The cumulative normal distribution function
Public Function CND2(x As Double) As Double
    
    If x = 0 Then
        CND2 = 0.5
    Else
        Dim L As Double, k As Double
        Const a1 = 0.31938153:  Const a2 = -0.356563782: Const a3 = 1.781477937:
        Const a4 = -1.821255978:  Const a5 = 1.330274429
    
        L = Abs(x)
        k = 1 / (1 + 0.2316419 * L)
        CND2 = 1 - 1 / Sqr(2 * Pi) * Exp(-L ^ 2 / 2) * (a1 * k + a2 * k ^ 2 + a3 * k ^ 3 + a4 * k ^ 4 + a5 * k ^ 5)
    
        If x < 0 Then
            CND2 = 1 - CND2
        End If
    End If
    
End Function


Public Function CBND4(A As Double, b As Double, rho As Double) As Double
    'modified/corrected from the second function in Drez & Wes paper pg. 105
    '0/0 case resolved by l'H rule

    Dim i As Integer
    Dim x As Variant, W As Variant
    Dim h1 As Double, h2 As Double
    Dim LH As Double, h12 As Double, h3 As Double, h5 As Double, h6 As Double, h7 As Double, h8 As Double
    Dim r1 As Double, r2 As Double, r3 As Double, rr As Double
    Dim AA As Double, ab As Double
  
    x = Array(0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992)
    W = Array(0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042)
    
    h1 = A
    h2 = b
    h12 = (h1 * h1 + h2 * h2) / 2
  
    If Abs(rho) >= 0.7 Then
        r2 = 1 - rho * rho
        r3 = Sqr(r2)
        If rho < 0 Then h2 = -h2
        h3 = h1 * h2
        h7 = Exp(-h3 / 2)
        If Abs(rho) < 1 Then
            h6 = Abs(h1 - h2)
            h5 = h6 * h6 / 2
            h6 = h6 / r3
            AA = 0.5 - h3 / 8
            ab = 3 - 2 * AA * h5
            LH = 0.13298076 * h6 * ab * (1 - CND(h6)) - Exp(-h5 / r2) * (ab + AA * r2) * 0.053051647
            For i = 0 To 4
                r1 = r3 * x(i)
                rr = r1 * r1
                r2 = Sqr(1 - rr)
                If h7 = 0 Then
                    h8 = 0
                Else
                    h8 = Exp(-h3 / (1 + r2)) / r2 / h7
                End If
                LH = LH - W(i) * Exp(-h5 / rr) * (h8 - 1 - AA * rr)
            Next i
        End If
        CBND4 = LH * r3 * h7 + CND(Min(h1, h2))
        If rho < 0 Then
            CBND4 = CND(h1) - CBND4
        End If
    Else
        h3 = h1 * h2
        If rho <> 0 Then
            For i = 0 To 4
                r1 = rho * x(i)
                r2 = 1 - r1 * r1
                LH = LH + W(i) * Exp((r1 * h3 - h12) / r2) / Sqr(r2)
            Next i
        End If
        CBND4 = CND(h1) * CND(h2) + rho * LH
    End If
     
End Function

Public Function CBNDGeneral(TypeFlag As Integer, x As Double, y As Double, rho As Double) As Double

    If TypeFlag = 1 Then 'Drezner-78
        CBNDGeneral = CBND2(x, y, rho)
    ElseIf TypeFlag = 2 Then 'Drezner-Weso-90a
        CBNDGeneral = CBND3(x, y, rho)
     ElseIf TypeFlag = 3 Then 'Drezner-Weso-90a
        CBNDGeneral = CBND4(x, y, rho)
     ElseIf TypeFlag = 4 Then ' Genze
        CBNDGeneral = CBND(x, y, rho)
    End If
    
End Function

' The  bivariate normal distribution function
Public Function BND(x As Double, y As Double, rho As Double) As Double
    
    BND = 1 / (2 * Pi * Sqr(1 - rho ^ 2)) * Exp(-1 / (2 * (1 - rho ^ 2)) * (x ^ 2 + y ^ 2 - 2 * x * y * rho))

End Function

Private Function ArcSin(x As Double) As Double
    If Abs(x) = 1 Then
        ArcSin = Sgn(x) * Pi / 2
    Else
        ArcSin = Atn(x / Sqr(1 - x ^ 2))
    End If
End Function


' The cumulative bivariate normal distribution function
' Based on Drezner-1978
Public Function CBND3(A As Double, b As Double, rho As Double) As Double

    Dim x As Variant, y As Variant
    Dim rho1 As Double, rho2 As Double, delta As Double
    Dim a1 As Double, b1 As Double, sum As Double
    Dim i As Integer, j As Integer
    
    x = Array(0.24840615, 0.39233107, 0.21141819, 0.03324666, 0.00082485334)
    y = Array(0.10024215, 0.48281397, 1.0609498, 1.7797294, 2.6697604)
    a1 = A / Sqr(2 * (1 - rho ^ 2))
    b1 = b / Sqr(2 * (1 - rho ^ 2))
    
    If A <= 0 And b <= 0 And rho <= 0 Then
        sum = 0
        For i = 0 To 4
            For j = 0 To 4
                sum = sum + x(i) * x(j) * Exp(a1 * (2 * y(i) - a1) + b1 * (2 * y(j) - b1) + 2 * rho * (y(i) - a1) * (y(j) - b1))
            Next
        Next
        CBND3 = Sqr(1 - rho ^ 2) / Pi * sum
    ElseIf A <= 0 And b >= 0 And rho >= 0 Then
        CBND3 = CND(A) - CBND3(A, -b, -rho)
    ElseIf A >= 0 And b <= 0 And rho >= 0 Then
        CBND3 = CND(b) - CBND3(-A, b, -rho)
    ElseIf A >= 0 And b >= 0 And rho <= 0 Then
        CBND3 = CND(A) + CND(b) - 1 + CBND3(-A, -b, rho)
    ElseIf A * b * rho > 0 Then
        rho1 = (rho * A - b) * Sgn(A) / Sqr(A ^ 2 - 2 * rho * A * b + b ^ 2)
        rho2 = (rho * b - A) * Sgn(b) / Sqr(A ^ 2 - 2 * rho * A * b + b ^ 2)
        delta = (1 - Sgn(A) * Sgn(b)) / 4
        CBND3 = CBND3(A, 0, rho1) + CBND3(b, 0, rho2) - delta
    End If
    
End Function
