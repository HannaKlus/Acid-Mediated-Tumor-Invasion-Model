#Phase 3: Clinical Simulation & Treatment Response
import numpy as np
import matplotlib.pyplot as plt

#-FUNCTION FOR SIMULATION-
def med_simulation (t_med, chem_s, duration=None):
    #Configuration
    dt = 0.01
    t_max = 250
    N = 300
    L = 300
    dx = L / N
    t = np.arange(0, t_max, dt)
    Nt = len(t)
    #
    u = np.zeros([Nt, N])
    v = np.zeros([Nt, N])
    w = np.zeros([Nt, N])
    history = np.zeros(Nt)
    #
    D_v = 0.5
    D_w = 20.0
    delta_pre = 1.2
    delta_post = 0.5
    a_decay = 0.5
    r_v = 0.8
    r_w = 1.0
    #
    u[0, :] = 1.0
    w[0, :] = 0.0
    mid = N // 2
    v[0, mid - 10:mid + 10] = 0.1
    #-NUMERICAL INTEGRATION (Euler Method)-
    for i in range(Nt - 1):
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

        U = u[i]
        V = v[i]
        W = w[i]

        lap_w = np.zeros(N)
        lap_v = np.zeros(N)
        lap_w[1:-1] = (W[:-2] - 2 * W[1:-1] + W[2:]) / dx ** 2
        lap_v[1:-1] = (V[:-2] - 2 * V[1:-1] + V[2:]) / dx ** 2

        # Reactions
        du = U * (1 - U) - delta * U * W
        dv = r_v * V * (1 - V) - chem * V
        dw = r_w * V - a_decay * W

        # Update
        u[i + 1] = U + du * dt
        v[i + 1] = V + (dv + D_v * lap_v) * dt
        w[i + 1] = W + (dw + D_w * lap_w) * dt

        # Safety
        u[i + 1, u[i + 1] < 0] = 0
        v[i + 1, v[i + 1] < 0] = 0

        history[i] = np.sum(v[i]) * dx

    history[-1] = np.sum(v[-1]) * dx
    return t, history

#-SCENARIOS-
# success
t1, h1 = med_simulation(t_med=40, chem_s=2.5)
# failure
t2, h2 = med_simulation(t_med=120, chem_s=0.5, duration=40)
# relapse
t3, h3 = med_simulation(t_med=120, chem_s=2.5, duration=40)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
fig.suptitle( 'Comparison of Clinical Outcomes', fontsize= 15)

axes = [ax1, ax2, ax3]
for ax in axes:
    ax.grid(axis='x', color='#e0e0e0', linestyle='--', linewidth=1)

# Scenario A (Relapse)
ax1.plot(t3, h3, color='#E03C39', linewidth=3, solid_capstyle='round')
ax1.text(t3[-1] + 2, h3[-1], "Relapse (Recurrence)", color='#E03C39', fontweight='bold', va='center')
ax1.axvspan(120, 160, color='#E03C39', alpha=0.1)

# Scenario B (Failure)
ax2.plot(t2, h2, color='#FF913A', linewidth=3, solid_capstyle='round')
ax2.text(t3[-1] + 2, h3[-1], "Failure", color='#FF913A', fontweight='bold', va='center')
ax2.axvspan(120, 160, color='#FF913A', alpha=0.1)
# Scenario C (Success)
ax3.plot(t1, h1, color='#DAA520', linewidth=3, solid_capstyle='round')
ax3.text(t1[-1] + 2, h1[-1], "Success", color='#DAA520', fontweight='bold', va='center')
ax3.axvspan(40, 300, color='#DAA520', alpha=0.1)


plt.tight_layout()
plt.subplots_adjust(hspace=0.05)
plt.savefig('glioma_results_med.png', dpi=300, bbox_inches='tight')
plt.show()





