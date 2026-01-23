# Acid-Mediated-Tumor-Invasion-Model
Acid-Mediated Tumor Invasion Mathematical Model Created in Python.
# Glioma growth and treatment simulation
A mathematical simulation of Glioma growth and treatment response, developed in Python.
This project is based on the seminal paper **"A Reaction-Diffusion Model of Cancer Invasion"** by Robert A. Gatenby and Edward T. Gawlinski.

The simulation was built in three evolutionary stages to progressively model the complexity of tumor dynamics.

##Project evolution
### Phase 1: Temporal Dynamics (ODE)
I started by modeling the population growth over time using the **Logistic Growth Equation**.
* *Goal:* Understand carrying capacity and basic growth limits.
*  *Mathematics:*Ordinary Differential Equations (ODE) computed using Euler's method.
![Phase 1 Result](img/glioma_results_ode.png)
*(Fig 1. Tumor cell population growth with medication showing full recovery)*
---
### Phase 2: Spatiotemporal Dynamics (PDE)
I introduced space into the model using the *Fisher-Kolmogorov Equation* (Reaction-Diffusion).
* *Goal:* Simulate how the tumor (and acidity) spreads into physical space
*  *Mathematics:* Partial Differential Equations (PDE) solved via Finite Difference Method.
![Phase 2 Result](img/glioma_results_pde.png)
*(Fig 2. Traveling wave showing tumor invasion into healthy tissue and tumor mass over time)*
---
### Phase 3: Clinical Simulation
The final model incorporates **chemotherapy, drug resistance, and relapse mechanisms**.
My goal was to demonstrate how changing two key parameters—**chemotherapy strength** and **duration of treatment**—drastically affects the clinical outcome.

![Phase 3 Result](img/glioma_results_med.png)
*(Fig 3. Comparison of Cure, Failure, and Relapse scenarios)*

## Tech Stack 
* *Python 3.x*
* *NumPy* (Numerical computing)
* *Matplotlib* (Data visualization)

## How to Run
1. Clone the repository.
2. Install dependencies:
   ```bash
   pip install numpy matplotlib
