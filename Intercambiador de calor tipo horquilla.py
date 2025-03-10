import numpy as np
import math as mt

#Diseño de intercambiador de tubos concéntricos tipo horquilla
#1cP = 2.4191 #lb/h*ft

#Datos del fluído caliente (Gasolina) Externo
Ti = 160 #°F
To = 120 #°F
Cph = 0.55 #Btu/lbm °F
kh = 0.087 #Btu/h ft °F
muh = 0.45 #cP
sh = 0.74
Prh = 6.88
Rdh = 0.001 #h ft2 °F/Btu

#Datos del fluído frío (Kerosene) Interno
mc = 15000 #lb/h
ti = 75 #°F
to = 120 #°F
Cpc = 0.48 #Btu/lbm °F
kc = 0.081 #Btu/h ft °F
muc = 1.5 #cP
sc = 0.82
Prc = 21.5
Rdc = 0.001 #h ft2 °F/Btu

#Cálculo de flujo másico fluído caliente
mh = mc*Cpc*(to-ti)/(Cph*(Ti-To))

#Temperatura media logarítmica
Tlm = ((Ti-to)-(To-ti))/(np.log((Ti-to)/(To-ti)))

#Parámetros del tubo
L = 16 #ft
k = 26 #Btu/h ft °F

#Tubo interno
#2 in SH40
di = 2.067/12 #ft
do = 2.375/12 #ft

#Tubo externo
#3 in SH40
Di = 3.068/12 #ft

#Cálculos térmicos
Rei = (4*mc)/(np.pi*di*muc*2.419)
if Rei < 2100:
  print('El flujo interno del intercambiador es de régimen laminar')
elif 10**4 > Rei > 2100:
  print('El flujo interno del intercambiador es de régimen transitorio')
else:
  print('El flujo interno del intercambiador es de régimen turbulento')
Afo = (np.pi/4)*(Di**2-do**2) #ft2

#Cálculo de flujo anular
De = (Di-do)
Go = mh/Afo #lb/h ft2
Reo = De*Go/(muh*2.419)
if Reo < 2100:
  print('El flujo anular del intercambiador es de régimen laminar')
elif 10**4 > Reo > 2100:
  print('El flujo anular del intercambiador es de régimen transitorio')
else: 
  print('El flujo anular del intercambiador es de régimen turbulento')

RDe = Di/do


#Cálculo de flujo
def convective_interal(Re, Pr, K, D, L):
    if Re < 2100:
        h = (K/D)*1.86*(Re*Pr*D/L)**(1/3)#Btu/h ft2 °F
    elif 10**4 > Re > 2100:
        h = (K/D)*0.116*(Re**(2/3)-125)*Pr**(1/3)*(1+(D/L)**(2/3)) #Btu/h ft2 °F
    else: 
        h = (K/D)*0.023*Re**0.8*Pr**(1/3) #Btu/h ft2 °F
    return h

def convective_annulus(Re, Pr, K, D, L, RDe):
    if Re < 2100:
        h = (K/D)*3.66+1.2*(RDe)**0.8+(0.19*(1+0.14*(RDe)**0.5)*(Re*Pr*D/L)**0.8)/(1+0.117*(Re*Pr*D/L)**0.467)#Btu/h ft2 °F
    elif 10**4 > Re > 2100:
        h = (K/D)*0.116*(Re**(2/3)-125)*Pr**(1/3)*(1+(D/L)**(2/3)) #Btu/h ft2 °F
    else:
        h = (K/D)*0.023*Re**0.8*Pr**(1/3) #Btu/h ft2 °F
    return h


hi = convective_interal(Rei, Prc, kc, di, L)

ho = convective_annulus(Reo, Prh, kh, De, L, RDe)


#Cálculo de coeficiente global
Ud = ((do/(di*hi))+((do*np.log(do/di))/(2*k))+(1/ho)+(do*Rdc/di)+Rdh)**-1
print(f'El coeficiente global del intercambiador es: {Ud:.2f} Btu/h*ft2*°F')


q = mc*Cpc*(to-ti)

A = q/(Ud*Tlm)

A1hp = np.pi*do*L*2

Nhp = A/A1hp
Nhp = mt.ceil(Nhp)
Lt = Nhp*L*2

#Pérdidas internas
fi = 0.36733*Rei**-0.2314
Afi = (np.pi/4)*(di**2) #ft2
Gi = mc/Afi #lb/h ft2
DPfi = (fi*Lt*Gi**2)/(7.5*10**12*di*sc)
Dpri = 1.6*10**-13*(2*Nhp-1)*(Gi**2/sc)


#Pérdidas internas totales
Dpti = DPfi+Dpri
#Pérdidas externas
fo = 0.36733*Reo**-0.2314
DPfo = (fo*Lt*Go**2)/(7.5*10**12*De*sh)
Dpro = 1.6*10**-13*(2*Nhp-1)*(Go**2/sh)
#Pérdidas internas totales
Dpto = DPfo+Dpro


#Pérdidas totales
Dptotal = Dpti+Dpto

print(f'Las caídas de presión son de: {Dptotal:.2f} psia')



def F_general(R, x, P):
  if R==1:
    numerator = P * (1 - x)
    denominator = x * (1 - P) * np.log(((1 - x) / ((1 - P) ** (1/x))) + x)
    return numerator / denominator
  else:
    numerator = (R - x) * np.log((1 - P) / (1 - P * R))
    denominator = x * (R-1) * np.log(((R - x) / (R * (1 - P * R) ** (1/x))) + (x / R))
    return numerator / denominator

P = (Ti-To)/(To-ti)
print(P)
R = (to-ti)/(Ti-To)
print(R)
x = 2
Y = R*P
print(Y)
#F = F_general(R, x, P)
#print(F)