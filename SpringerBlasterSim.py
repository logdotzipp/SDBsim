import numpy as np
import math
import matplotlib.pyplot as plt

# Define variables
T_p  = 125/1000 #Plunger Travel (m)
L_b = 400/1000 #Barrel Length (m)
A_p  = (np.pi*((1.375*25.4)/2)**2)/(1000**2) #Cross Sectional Area of the Plunger Tube (m^2)
m_p = 100/1000 # Mass of the plunger (kg)
s = 10/1000 # Spring Preload (m) - Preload distance on the spring at rest (before the plunger is primed)
             # The total initial spring displacement is calculated as the Spring Preload + Plunger travel
K_s = 3.5 * 175.127 #Spring Rate (N/m)
A_b  = (np.pi*(12.7/2)**2)/(1000**2) #Cross Sectional Area of the barrel (m^2)
m_d = 1/1000 # Mass of the dart (kg)
V_d1 = 0.0000001 # Dead space volume on plunger side of orifice (m^3)
V_d2 = ((np.pi*(10/2)**2)*100)/(1000**3) # Dead space volume on barrel side of orifice (m^3)
Patm = 101.3 * 1000 # Atmospheric pressure (Pa aka N/m^2)
Tatm = 21 + 273.3 # Temp of air in plunger (K). Assumed constant.


# Orifice Discharge Flow coefficient. This is based on the geometry of the orifice (square vs fillet edge)
# Perfectly square edge corners result in a discharge coeff of around 0.5
# Significant Rounding of the orifice edge can increase the discharge coeff to about 0.94
# See Figure 10 from https://ntrs.nasa.gov/api/citations/19690028630/downloads/19690028630.pdf
C_d = 0.5
A_o = (np.pi*(10/2)**2)/(1000**2) # Area of the orifice (Smallest area in flow)

R = 287 # Ideal Gas constant for air (J/(kg*K)
k = 1.4 # Specific heat ratio of air

# The plunger friction and dart friction are modeled as constant forces. There is some normal force N that
# is applied in a circle outwards of the dart/plunger. This force is unknown, and similarly so is the
# coefficient of friction. It's also likely that increased pressure in the plunger tube expands the dart
# and oring, increasing the friction. The plunger friction is likely fairly negligible. However, the dart
# friction is somewhat critical to allowing pressure to build in the plunger tube.
# Further testing must be done to find an emprirical value for the frictional force of the dart in barrel.

# This model has static and dynamic friction as the same.
F_fp = 0.25 # Plunger Frictional Force (N)
F_fd = 2 # Dart/barrel Frictional Force (N)


# Calculate the initial mass of the air in each cylinder using ideal gas law
m1_0 = Patm*(A_p*T_p + V_d1)/(R*Tatm)
m2_0 = Patm*(V_d2)/(R*Tatm)

def mdot_calc (p1,p2):
    # Calculate the mass flowrate from the piston to the barrel
    # The size of the orifice limits the flowrate between the two cylinders
    # Use the isentropic model for compressible subsonic flow through an orifice

    a = (2/(k-1))*((p2/p1)**(2/k) - (p2/p1)**((k+1)/k))

    if a < 0:
        mdot = 0
    else:
        mdot = C_d*A_o*p1*math.sqrt(k/(R*Tatm))*math.sqrt(a)
    #print(mdot)
    #print(" ")
    return mdot

def p1_calc (m1, x1):
    # Calculate the pressure inside cylinder 1 using the ideal gas law
    V1 = (T_p - x1)*A_p + V_d1
    if V1 > 0:
        p1 = m1 * R * Tatm / V1
    else:
        p1 = Patm
    return p1

def p2_calc (m2, x2):
    # Calculate the pressure inside cylinder 2 using the ideal gas law
    V2 = A_b*x2 + V_d2
    if V2 > 0:
        p2 = m2 * R * Tatm / V2
    else:
        p2 = Patm
    return p2

def Function(t, u):
    x1, vx1, x2, vx2, m1, m2 = u  # unpack state vector
    # Calculate the pressure in each cylinder
    p1 = p1_calc(m1, x1)
    p2 = p2_calc (m2, x2)

    # Calculate the mass flow rate of air between the two cylinders
    mdot = mdot_calc(p1,p2)

    dm1dt = -mdot # Flow is out of cylinder 1
    dm2dt = mdot # Flow is into cylinder 2

    dx1dt = vx1
    dvx1dt = ((T_p + s)*K_s - K_s*x1 - (p1 - Patm)*A_p - F_fp)/m_p
    dx2dt = vx2
    dvx2dt = ((p2 - Patm)*A_b - F_fd)/m_d
    if(dvx2dt < 0): # Account for static friction
        dvx2dt = 0 

    return np.array([dx1dt, dvx1dt, dx2dt, dvx2dt, dm1dt, dm2dt])

def rk4_step(F, t, u, h):
    k1 = F(t, u)
    k2 = F(t + h/2, u + h*k1/2)
    k3 = F(t + h/2, u + h*k2/2)
    k4 = F(t + h, u + h*k3)
    return u + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

def rk4_solver(Func, u0, p0, t0, t_end, h):
    t_range = np.arange(t0, t_end + h, h)
    u_values = [u0]
    p_values = [p0]
    u = u0
    p = p0
    t_values = [0]
    

    for t in t_range[:-1]:
        u = rk4_step(Func, t, u, h)
        u_values.append(u)
        t_values.append(t)
        
        # Calculate pressures
        x1, vx1, x2, vx2, m1, m2 = u  # unpack state vector
        p1 = p1_calc(m1, x1) - Patm
        p2 = p2_calc (m2, x2) - Patm
        p_values.append(np.array([p1,p2]))

        if(x1 >= (T_p) or m1 <= 0): # Bottom out of plunger, end sim
            break
        
        if(x2 >= L_b): # Dart has exited barrel, end sim
            break
    
    return np.array(t_values), np.array(u_values), np.array(p_values)




u0 = np.array([0.0, 0.0, 0.0, 0.0, m1_0, m2_0])  # x1, vx1, x2, vx2, m1, m2 initial conditions
p0 = np.array([0.0, 0.0])
t, u, p = rk4_solver(Function, u0, p0, t0=0.0, t_end=0.1, h=0.00001)

x1, vx1, x2, vx2, m1, m2 = u.T  # unpack results
p1, p2 = p.T # unpack pressure data

x1_mm = x1 * 1000
vx1_mm = vx1 * 1000
x2_mm = x2 * 1000
vx2_mm = vx2 * 1000

vx2_ftpersec = vx2_mm * (1/(25.4*12))

p1_psi = p1 * 0.000145038
p2_psi = p2 * 0.000145038

print(p1)

print("")

if(x1[-1] >= T_p):
    # Plunger bottom out condition. Dart has not escaped barrel
    print("Dart Displacement (mm): " + str(x2_mm[-1]))
else:
    # Dart escaped barrel condition
    print("Exit FPS: " + str(vx2_ftpersec[-1]))
    x1_rem = (T_p - x1[-1])*1000
    print("Remianing plunger travel (mm): " + str(x1_rem))
print("")


# Create a figure with multiple subplots
fig, ax = plt.subplots(5, 1, figsize=(8, 8), sharex=True)

# Plot position over time
ax[0].plot(t, x1_mm, label='Plunger Displacement x1(t)', color='blue')
ax[0].plot(t, x2_mm, label='Dart Displacement x2(t)', color='orange', linestyle='--')
ax[0].set_ylabel('Position (mm)')
ax[0].legend()
ax[0].grid(True)

# Plot velocity over time
ax[1].plot(t, vx1_mm, label='Plunger Velocity vx1(t)', color='purple')
ax[1].plot(t, vx2_mm, label='Dart Velocity vx2(t)', color='green', linestyle='--')
ax[1].set_xlabel('Time [s]')
ax[1].set_ylabel('Velocity (mm/s)')
ax[1].legend()
ax[1].grid(True)

ax[2].plot(t, vx2_ftpersec, label='Dart Velocity vx2(t)', color='orange')
ax[2].set_xlabel('Time [s]')
ax[2].set_ylabel('Velocity (ft/s)')
ax[2].legend()
ax[2].grid(True)

ax[3].plot(t, p1_psi, label='Plunger psi1(t)', color='green')
ax[3].plot(t, p2_psi, label='Plunger psi2(t)', color='blue')
ax[3].set_xlabel('Time [s]')
ax[3].set_ylabel('Pressure [psi]')
ax[3].legend()
ax[3].grid(True)

ax[4].plot(t, m1, label='Plunger air m1(t)', color='red')
ax[4].plot(t, m2, label='Barrel air m2(t)', color='green')
#ax[4].plot(t, m1+m2, label='mtotal', color='orange')
ax[4].set_xlabel('Time [s]')
ax[4].set_ylabel('Mass (kg)')
ax[4].legend()
ax[4].grid(True)

plt.tight_layout()
plt.show()