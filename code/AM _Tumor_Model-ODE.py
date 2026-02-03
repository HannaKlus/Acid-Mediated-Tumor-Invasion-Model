#Acid-Mediated Tumor Invasion Model
#Phase: 1 - Temporal Dynamics (ODE Model)
# Modeling the interaction between Host Tissue, Tumor, and Acid.
import numpy as np
import matplotlib.pyplot as plt

# -TIME CONFIGURATION-
dt = 0.05
t_max = 100
t = np.arange(0, t_max, dt)

# -INITIALIZATION-
u = np.zeros(len(t)) #Host Tissue
v = np.zeros(len(t)) #Tumor Tissue
w = np.zeros(len(t)) #Acid Concentration

#-MODEL PARAMETERS-
#Interaction coefficients
delta_pre = 0.8 #Acid toxicity rate
delta_post = 0.4
a_decay = 0.5 #Natural acid decay/removal rate
#Growth rates
r_v = 0.8 #Tumor proliferation rate
r_w = 1.0 #Acid production rate
#Treatment parameters
t_med = 10 #Time when we start with medication
chem_s = 1.0 #Chemotherapy strength
duration = 80

#-INITIAL CONDITIONS (t=0)
u[0] = 1.0
v[0] = 0.01
w[0] = 0.0

#-NUMERICAL INTEGRATION (Euler Method)-
for i in range(0, len(t)-1):
    #Treatment switch logic
    current_time = t[i]
    chem_active = False

    if duration is None:
        if current_time >= t_med:
            chem_active = True
    else:
        if current_time >= t_med and current_time < (t_med + duration):
            chem_active = True

    if chem_active:
        delta = delta_post
        chem = chem_s
    else:
        delta = delta_pre
        chem = 0.0

    #Current states
    U = u[i]
    V = v[i]
    W = w[i]
    #Differential equations
    du = U*(1-U) - delta * U * W
    dv = r_v * V*(1-V) - chem * V
    dw = r_w * V - a_decay * W
    #Euler Step
    u[i+1] = u[i] + dt * du
    v[i+1] = v[i] + dt * dv
    w[i+1] = w[i] + dt * dw
    #Extinction Thresholds (Prevention of "Zombie" effect)
    if u[i+1] < 0.02:
        u[i+1] = 0.0
    #Non-negative constraint
    if v[i+1] < 0.0:
        v[i+1] = 0.0
    if w[i+1] < 0.0:
        w[i+1] = 0.0

#-PLOTTING-
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Population Dynamics
ax1.plot(t, u, label = 'Healthy Tissue', linewidth = 2)
ax1.plot(t, v, label = 'Tumor Tissue', linewidth = 2)

ax1.set_xlabel('Time [arbitrary units]')
ax1.set_ylabel('Normalized Density')
ax1.set_title('Population Dynamics over Time')
#Treatment highlighting
ax1.axvline(x = t_med, color = 'red', linestyle = '--', linewidth = 1, alpha = 0.5)
ax1.axvline(x = t_med+duration, color = 'red', linestyle = '--', linewidth = 1, alpha = 0.5)
ax1.axvspan(t_med, t_med+duration, color = 'red', alpha = 0.1)
ax1.legend()

# Plot 2: Phase Portrait
ax2.plot(v, u, linewidth = 2, color = 'purple')
ax2.set_title('Phase Portrait: Tumor-Host Interaction')
ax2.set_ylabel('Healthy Tissue Density')
ax2.set_xlabel('Tumor Tissue Density')
plt.grid(True)
plt.tight_layout()
plt.savefig('glioma_results_ode.png', dpi=300, bbox_inches='tight')
plt.show()




