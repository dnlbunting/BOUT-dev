#
# Generate the test case using SymPy
#
# Equations are:
#
# 
# Delp2(phi) = vort + Sphi

from boutdata.mms import *

####

#from sympy import Integral
#def rms(expr):
#    """ Calculates the RMS value of a function over space
#    """
#    ix = Integral(expr**2, (x, 0, 1))
#    ixy = Integral(ix, (y, 0, 2*pi))
#    ixyz = 
#    return sqrt( , (y, 0, 2*pi), (z, 0, 2*pi)).evalf() )
    

def C(f):
    """ Curvature 
    """
    return bxcvz * diff(f, metric.z)

#######################################
# Inputs, essentially as in BOUT.inp

# Switches
ionvis  = False # Ion Viscosity
Ti      = 10    # Ion temperature for viscosity calculation
elecvis = False # Electron viscosity
resistivity = True

parallel = False  # Parallel dynamics

# Parameters
Tnorm = 5     # Electron Temperature (eV)
Nnorm = 2e18  # Background plasma density (m^-3)
Bnorm = 0.35  # Magnetic field [T]
AA    = 2     # Ion atomic mass

#######################################
# 

mi_me = AA * 1.67262158e-27 / 9.109e-31
beta_e = qe*Tnorm*Nnorm / (Bnorm**2/mu0)

# Normalisation parameters
Cs0      = sqrt(qe*Tnorm / (AA*Mp))
Omega_ci = qe*Bnorm / (AA*Mp)
rho_s0   = Cs0 / Omega_ci

# Define a metric

Rxy = 1.5    # Major radius
Bpxy = 0.35  # Poloidal magnetic field
Bxy = 0.35   # Total magnetic field
Btxy = 0.0   # Toroidal magnetic field
hthe = 1.0   # Poloidal arc length

sinty = 0.0
sbp = 1.0

dx = 2e-5

nx = 130
ZMAX = 1e-3

MXG = 1
MYG = 0

bxcvz = 100 # Curvature

# Normalise

Rxy /= rho_s0
Bpxy /= Bnorm
Btxy /= Bnorm
Bxy  /= Bnorm
hthe /= rho_s0

dx   /= rho_s0**2 * Bnorm

bxcvz *= rho_s0**2;

# Define a metric
metric = Metric() # Identity

Lx = dx * (nx - 2.*MXG) # Size of the X domain

metric.g11 = (Rxy*Bpxy)**2
metric.g22 = 1.0 / (hthe**2)
metric.g33 = (sinty**2)*metric.g11 + (Bxy**2)/metric.g11
metric.g12 = 0.0
metric.g13 = -sinty*metric.g11
metric.g23 = -sbp*Btxy/(hthe*Bpxy*Rxy)
  
metric.J = hthe / Bpxy
B = metric.B = Bxy
  
metric.g_11 = 1.0/metric.g11 + ((sinty*Rxy)**2)
metric.g_22 = (Bxy*hthe/Bpxy)**2
metric.g_33 = Rxy*Rxy
metric.g_12 = sbp*Btxy*hthe*sinty*Rxy/Bpxy
metric.g_13 = sinty*Rxy*Rxy
metric.g_23 = sbp*Btxy*hthe*Rxy/Bpxy


# Define the manufactured solution in terms of input x,y and z

Ne    = 0.9 + 0.9*x + 0.2*cos(t)*sin(5.*x**2 - 5.*z)
Te    = 1 + 0.5*cos(t)*cos(3*x**2 - 4*z)
Vort  = sin(2*t)*cos(x - z)
VePsi = 0.5*cos(7*t)*cos(3*x**2 - 3*z)
Vi    = -0.1*cos(7*t)*cos(3*x**2 - 3*z)

phi = sin(z - x + t)*sin(2.*pi*x) # Must satisfy Dirichlet BCs for now
psi = sin(pi*x) *(0.5*x - cos(7*t)*sin(3.*x**2 - 3*z)) # Must satisfy Dirichlet BCs for now

# Substitute to get in terms of actual x,y,z coordinates

replace = [ (x, metric.x / Lx), (z, metric.z / ZMAX) ]

Ne    = Ne.subs(replace)
Te    = Te.subs(replace)
Vort  = Vort.subs(replace)
VePsi = VePsi.subs(replace)
Vi    = Vi.subs(replace)
phi   = phi.subs(replace)
psi   = psi.subs(replace)

# Calculate RHS function

Pe = Ne * Te
Ve = VePsi - 0.5*mi_me*beta_e*psi

Gi = 0.0
Ge = 0.0

nu = 0.0

### Electron density

dNedt = (
    - bracket(phi, Ne, metric)
    + (2/B) * ( C(Pe) - Ne*C(phi) )
    )

if parallel:
    dNedt -= Ne*Grad_par(Ve, metric) + Vpar_Grad_par(Ve, Ne, metric)

### Electron temperature

dTedt = (
    - bracket(phi, Te, metric)
    + (4./3)*(Te/B) * ( (7./2)*C(Te) + (Te/Ne)*C(Ne) - C(phi) )
    )

if parallel:
    dTedt -= Vpar_Grad_par(Ve, Te, metric);
    dTedt += (2./3.)*Te*( 0.71*Grad_par(Vi, metric) - 1.71*Grad_par(Ve, metric) + 0.71*(Vi-Ve)*Grad_par(log(Ne), metric))

### Vorticity

dVortdt = (
    - bracket(phi, Vort, metric)
    + 2.*B*C(Pe)/Ne + B*C(Gi)/(3.*Ne)
    )

if parallel:
    dVortdt += (
        - Vpar_Grad_par(Vi, Vort, metric)
        + B^2 * (Grad_par(Vi - Ve, metric) + (Vi - Ve)*Grad_par(log(Ne), metric))
        )

### Parallel Ohm's law

if parallel:
    dVePsidt = (
        - bracket(phi, Ve, metric)
        - Vpar_Grad_par(Ve, Ve, metric)
        - mi_me*(2./3.)*Grad_par(Ge, metric)
        - mi_me*nu*(Ve - Vi)
        + mi_me*Grad_par(phi, metric)
        - mi_me*(  Te*Grad_par(log(Ne), metric) + 1.71*Grad_par(Te, metric) )
        )
else:
    dVePsidt = 0.0

### Parallel ion velocity

if parallel:
    dVidt = (
        - bracket(phi, Vi, metric)
        - Vpar_Grad_par(Vi, Vi, metric)
        - (2./3.)*Grad_par(Gi, metric)
        - (Grad_par(Te, metric) + Te*Grad_par(log(Ne), metric))
        )
else:
    dVidt = 0.0

#######################

# Create list of all evolving variables
vars = [
    (Ne, dNedt, "Ne"),
    (Te, dTedt, "Te"),
    (Vort, dVortdt, "Vort"),
    (VePsi, dVePsidt, "VePsi"),
    (Vi, dVidt, "Vi")
    ]

replace = [ (metric.x, x * Lx), (metric.z, z * ZMAX) ]

# Loop over variables and print solution, source etc.
for v, dvdt, name in vars:
    # Calculate source
    S = diff(v, t) - dvdt
    
    dvdx = diff(v, metric.x)
    dvdy = diff(v, metric.y)

    # Substitute back to get in terms of x,y,z
    v = v.subs(replace)
    dvdx = dvdx.subs(replace)
    dvdy = dvdy.subs(replace)
    S = S.subs(replace)
    
    print("["+name+"]")
    print("solution = " + exprToStr(v))
    print("\nddx     = " + exprToStr(dvdx))
    print("\nddy     = " + exprToStr(dvdy))
    print("\nsource   = " + exprToStr(S))

# Potential
Sphi = Delp2(phi, metric) - Vort

print("\n[phi]")
print("solution = " + exprToStr(phi.subs(replace)))
print("\nsource   = " + exprToStr(Sphi.subs(replace)))


print "\n\nDelp2 phi = ", Delp2(phi, metric).subs(replace)

##########################################
# Check magnitudes of individual terms

print "\n\nDensity terms:"
print "  bracket : ", (-bracket(phi, Ne, metric)).subs(replace), "\n"
print "  (2/B)*C(Pe) : ", ((2/B) * C(Pe)).subs(replace), "\n"
print "  (2/B)*Ne*C(phi) : ", ((2/B)*Ne*C(phi)).subs(replace), "\n"


print "\n\nTemperature terms:"
print "  bracket : ", (- bracket(phi, Te, metric)).subs(replace)
print "  C(Te)   : ", ((4./3)*(Te/B) *(7./2)*C(Te)).subs(replace)
print "  (Te/Ne)*C(Ne) : ", ((4./3)*(Te/B) *(Te/Ne)*C(Ne)).subs(replace)
print "  C(phi)  : ", (-(4./3)*(Te/B) * C(phi)).subs(replace)



print "\n\nVorticity terms:"
print "  bracket : ", (- bracket(phi, Vort, metric)).subs(replace)
print "  C(Pe)   : ", (2.*B*C(Pe)/Ne).subs(replace)
print "  C(Gi)   : ", B*C(Gi)/(3.*Ne)
