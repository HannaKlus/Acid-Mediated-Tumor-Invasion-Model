#Acid-Mediated Tumor Invasion Model
#Phase 2 - PDE Model

#-LIBRARIES-
import numpy as np
import matplotlib.pyplot as plt

#-TIME AND SPACE CONFIGURATION-
dt = 0.01
t_max = 200
t = np.arange(0, t_max, dt)
N = 200
L = 400
dx = L/N
x = np.linspace(0, L, N)
Nt = len(t)

#-INITIALIZATION-
u = np.zeros([Nt, N])
v = np.zeros([Nt, N])
w = np.zeros([Nt, N])
tumor_weight_history = np.zeros(Nt)

#-MODEL PARAMETERS-
#Diffusion constants
D_u = 0
D_v = 0.1
D_w = 40.0
#Other parameters
delta_pre = 3.0 #Acid toxicity rate
delta_post = 0.5
a_decay = 0.3 #Natural acid decay/removal rate
#Growth rates
r_v = 1.0 #Tumor proliferation rate
r_w = 1.1 #Acid production rate
#Treatment parameters
t_med = 60 #Time when we start with medication
chem_s = 0.5 #Chemotherapy strength

#-INITIAL CONDITIONS (t=0)
u[0, :] = 1.0
w[0, :] = 0.0
#
middle = N // 2
width = 10
v[0, middle-width : middle+width] = 0.1

#-NUMERICAL INTEGRATION (Euler Method)-
for i in range(Nt-1):
    if t[i]<t_med:
        delta = delta_pre
        chem = 0.0
    else:
        delta = delta_post
        chem = chem_s
    # Current states
    U = u[i, :]
    V = v[i, :]
    W = w[i, :]
    #Acid Diffusion
    lap_w = np.zeros(N)
    lap_w[1:-1] = (W[:-2] - 2 * W[1:-1] + W[2:]) / dx ** 2
    #Tumor Diffusion
    lap_v = np.zeros(N)
    lap_v[1:-1] = (V[:-2] - 2 * V[1:-1] + V[2:]) / dx ** 2
    #
    du = U * (1 - U) - delta * U * W
    dv = r_v * V * (1 - V) - chem * V
    dw = r_w * V - a_decay * W
    #
    u[i+1,:] = U + du * dt
    v[i+1,:] = V + (dv + lap_v * D_v) * dt
    w[i+1,:] = W + (dw + lap_w * D_w) * dt
    #
    u[i + 1, u[i + 1, :] < 0.3] = 0.0
    #
    current_weight = np.sum(v[i, :]) * dx
    tumor_weight_history[i] = current_weight

#-PLOTTING-
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
ax1.plot(x, u[-1, :], label='Healthy Tissue Density', color='blue')
ax1.plot(x, v[-1, :], label='Tumor Tissue Density', color='red')
ax1.plot(x, w[-1, :], label='Acid', color='green', linestyle='--')
ax1.set_title(f'Tissue state at T={t_max}')
ax1.set_xlabel('Space [a.u.]')
ax1.set_ylabel('Tissue Density')
ax1.fill_between(x, u[-1, :], color='blue', alpha=0.1)
ax1.fill_between(x, v[-1, :], color='red', alpha=0.1)
ax1.legend()

ax2.plot(t, tumor_weight_history, label='Total Tumor Weight', color='purple', linewidth=2)
ax2.axvline(x=t_med, color='black', linestyle=':', label='Start of medication')
ax2.fill_between(t, 0, tumor_weight_history, color='purple', alpha=0.1)
ax2.set_title('Tumor Growth History (Time)')
ax2.set_xlabel('Time [a.u.]')
ax2.set_ylabel('Total Mass')
ax2.legend()

plt.tight_layout()
plt.grid(True, alpha=0.3)
plt.savefig('glioma_results_pde.png', dpi=300, bbox_inches='tight')
plt.show()





