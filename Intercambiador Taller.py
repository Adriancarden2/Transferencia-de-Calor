import numpy as np
import sys

#Textos replicables
Re = f"Número de Reynolds es"
Cond = "El flujo interno del intercambiador es de régimen"
Cond2 = "El flujo externo del intercambiador es de régimen"

#PROPIEDADES DE FLUIDOS 
#Datos sobre acetato de amilo (Interno)
mc = 18500/2 #lbm/h
ti = 175 #°F
to = 125 #°F
Cpc = 0.52 #Btu/lbm°F
kc = 0.0615 #Btu/hft°F
muc = 1.1130 #lbm/ft h
sc = 0.83
Prc = Cpc*muc/kc
Rdc = 0.001 #hft2°F/Btu

#Datos sobre alcohol amílico (Externo)
Ti = 75 #°F
To = 120 #°F
Cph = 0.59 #Btu/lbm°F
kh = 0.0865 #Btu/hft°F
muh = 6.0730 #lbm/ft h
sh = 0.81
Prh = Cph*muh/kh
Rdh = 0.001 #hft2°F/Btu

#Datos del tubo
L = 20 #ft
k = 26 #Btu/hft°F
#Diametro interno (2 in)
di = 2.067/12 #ft
do = 2.375/12 #ft
#Diametro externo (3 in)
Di = 3.068/12 #ft

#Cálculos térmicos
mh = mc*(Cpc*(ti-to))/(Cph*(To-Ti))
R1 = mh*Cph*(To-Ti) - mc*Cpc*(ti-to)

if R1 < 0:
    sys.exit("ERROR en mh")
else: 
    R1==0
    print(mh)

#Temperatura media logarítmica
Tlm = ((ti-To)-(to-Ti))/np.log((ti-To)/(to-Ti))
print(Tlm)

#Reynolds interno
Rei = 4*mc/(np.pi*di*muc)
if Rei < 2100:
  print(Cond, 'laminar')
elif 10**4 > Rei > 2100:
  print(Cond, 'transitorio')
else:
  print(Cond, 'turbulento')

#Reynolds externo
Afo = (np.pi/4)*(Di**2-do**2) #ft2
De = (Di - do) #ft
Go = mh/Afo #lb/h.ft2
Reo = De*Go/(muh)
if Reo < 2100:
  print(Cond2, 'laminar')
elif 10**4 > Reo > 2100:
  print(Cond2, 'transitorio')
else:
  print(Cond2, 'turbulento')

print(Re, Rei)
print(Re, Reo)

#Cálculo de flujo
RDe = Di/do
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

#Suposición de longitud
Ls = 240 #ft

#Coeficientes convectivos
hi = convective_interal(Rei, Prc, kc, di, Ls)
print(f'El coeficiente convectivo interno es de {hi:.2f}')
ho = convective_annulus(Reo, Prh, kh, De, Ls, RDe)
print(f'El coeficiente convectivo externo es de {ho:.2f}')

#Cálculo de coeficiente global
Ud = ((do/(di*hi))+do*np.log(do/di)/(2*k)+(1/ho)+ (do*Rdc/di)+Rdh)**-1 #btu/h.ft2.F
print(f'El coeficiente global del intercambiador es: {Ud:.2f} Btu/h*ft2*°F')

DPtis = 1
DPtos = 100
#Factor de corrección
def P_correction(Ti, To, ti, to):
    if DPtis > DPtos:
       P = (to-ti)/(Ti-Ti)
    elif DPtis == DPtos:
       sys.exit("No genera cambios significativos hacer división de flujo")
    else:
       P = (To-Ti)/(ti-Ti)  
    return P

def R_correction(Ti, To, ti, to):
    if DPtis > DPtos:
       R = (ti-to)/(To-Ti)
    elif DPtis == DPtos:
       sys.exit("No genera cambios significativos hacer división de flujo")
    else:
       R = (Ti-To)/(to-Ti)  
    return R

P = P_correction(Ti, To, ti, to)
R = R_correction(Ti, To, ti, to)

def F_general(R, x, P):
  if R==1:
    numerator = P * (1 - x)
    denominator = x * (1 - P) * np.log(((1 - x) / ((1 - P) ** (1/x))) + x)
    return numerator / denominator
  else:
    numerator = (R - x) / (x * (R - 1))
    log1 = np.log((1 - P) / (1 - P * R))
    denominator = np.log(((R - x) / (R * (1 - P * R) ** (1/x))) + (x / R))
    return (numerator * log1) / denominator

x=2
F = F_general(R, x, P)
print(F)

q = mc*Cpc*(ti-to) # Btu/h
A = q/(Ud*Tlm*F) # Área total de intercambio de calor ft2
print(f'El área total de intercambio de calor es: {A:.2f} ft2')

# Área de una horquilla
A1hp = np.pi*do*L*2
#Número de horquillas
Nhp = A/A1hp
Nhp = np.ceil(Nhp)
Nhp = 6
print(f'El número de horquillas es: {Nhp}')
Lt = Nhp*L*2 #ft
print(f'La longitud total de las horquillas es: {Lt:.2f} ft')
#Cálculos hidráulicos

#Pérdias internas
fi = 0.36733*(Rei)**-0.2314
Afi = (np.pi/4)*(di**2) #ft2
Gi = mc/(Afi)
DPfi = fi*Lt*Gi**2/(7.5*10**12*di*sc) #psi
DPri = 1.6*10**-13*(2*Nhp-1)*Gi**2/sc #psi

#pédidas internas totales
DPti = DPfi + DPri #psi
print(f'Las pérdidas internas son de: {DPti:.2f} psi')

#pérdidas externas
fo = 0.36733*(Reo)**-0.2314
DPfo = fo*Lt*Go**2/(7.5*10**12*De*sh) #psi
DPro = 1.6*10**-13*(2*Nhp-1)*Go**2/sh #psi

#Pérdidas por boquillas
# 2 in ced 40
dn = 2.067/12 # ft
Afn = (np.pi/4)*(dn**2) #ft2
Gn = mh/(Afn)
DPfn = 2.0*10**-13*(Nhp*Gn**2/sh) #psi

#pédidas externas totales
DPto = DPfo + DPro + DPfn #psi
print(f'Las pérdidas externas son de: {DPto:.2f} psi')

#Pérdidas totales
DPtotal = DPti + DPto #psi
print(f'Las pérdidas totales son de: {DPtotal:.2f} psi')

#Cálculos económicos
Npipe = Lt/20 #Número de tubos
Npipe = np.ceil(Npipe)
print(f'El número de tubos es: {Npipe}')
cost2in = 52598 #COP
cost3in = 103887 #COP

cost_total = Npipe*cost2in + Npipe*cost3in
print(f'El costo total de los tubos es: {cost_total:.2f} COP')