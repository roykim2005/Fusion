import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar

config_data = """{"qc": 1.5, "qa": 3.6, "pc": 1.2e-1, "epsa": 3, "B0": 5.3, "Nr": 1000}"""
params = json.loads(config_data)
qa, pc, B0, Nr, epsa = params["qa"], params["pc"], params["B0"], params["Nr"], params["epsa"]
qc = 0.7

# Integrates d\psi^*/dr and finds \psi^*(r) = \int_0^r (1 - q(r')) B_\theta(r') dr'

mu = qa / qc
r = np.linspace(1e-5, 1.0, Nr) # r/a, normalized
q = qa * (r**2) / (1 - (1 - r**2)**mu) # safety factor from Wesson. was gonna use mu=2 but root finding algo failed
p = pc * (1.0 - r**2)**mu
B_theta = (r * B0 * epsa) / q

dpsi_star_dr = (1.0 - q) * B_theta # helical magnetic field

# integrates
psi_star = np.zeros_like(r)
for i in range(1, len(r)):
    psi_star[i] = np.trapezoid(dpsi_star_dr[:i+1], r[:i+1])

# q(r1) = 1 res surface, interpolates and finds r1
q_func = interp1d(r, q - 1.0, kind='cubic')
r1 = float(root_scalar(q_func, bracket=[r[0], r[-1]]).root)

psi_star_interp = interp1d(r, psi_star, kind='cubic')
psi_star_0 = psi_star[0]

# outer boundary where reconnection stops, psi_star(r0) = psi_star(0)
r_outer = r[r > r1]
psi_outer = psi_star[r > r1] - psi_star_0
r0_func = interp1d(psi_outer, r_outer, kind='cubic')
r0 = float(r0_func(0.0))

rminus = np.linspace(1e-5, r1, 300)
psi_1 = psi_star_interp(rminus)

r_domain = r[r >= r1]
psi_domain = psi_star[r >= r1]

psi_1_clamped = np.clip(psi_1, psi_domain.min(), psi_domain.max())

rplus_func = interp1d(psi_domain, r_domain, kind='cubic', fill_value="extrapolate")
rplus = rplus_func(psi_1_clamped)
r_0 = np.sqrt(rplus**2 - r1**2)

q_pre = q
q_r0 = float(interp1d(r, q_pre, kind='cubic')(r0))
q_post = np.copy(q_pre)
maskcore = r <= r0
q_post[maskcore] = 1.0 + (q_r0 - 1.0) * (r[maskcore]/(r0))**2 # tried to use np.where(r <= r0, 1.0, q_pre), and it has an annoying discontinuity problem at r0. 
# q_post satisfies q(0) = 1, q(r0) = q_r0, and dq/dr > 0. trying to see if there are better models. 


plt.figure(figsize=(8, 5))
plt.plot(r, psi_star, 'b-', label='Ψ*(r)')
plt.axvline(r1, color='r', linestyle='--', label=f'r_1={r1:.2f}')
plt.axvline(r0, color='g', linestyle='--', label=f'r_0={r0:.2f}')
plt.xlabel('r/a', fontsize = 16)
plt.ylabel('Ψ*(r)', fontsize = 16)
plt.title('r vs Ψ*(r)', fontsize = 16)
plt.grid(True)
plt.legend()
plt.show()

plt.plot(r, q_pre, 'k--', linewidth=2, label='Pre crash q(r)')
plt.plot(r, q_post, 'r-', linewidth=2.5, label='Post crash q(r)')
plt.axvline(r1, color='blue', linestyle=':', label=f'r_1 = {r1:.2f}')
plt.axvline(r0, color='green', linestyle=':', label=f'r_0 = {r0:.2f}')
plt.axhline(1.0, color='gray', linestyle=':', alpha=0.6)
plt.xlabel('r/a')
plt.ylabel('q(r)')
plt.title('Pre crash and post crash q(r) profiles')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

q_pre_func = interp1d(r, q_pre, kind='cubic')
q_post_func = interp1d(r, q_post, kind='cubic')

print(f"qc : {qc:.2f}")
print(f"r1 : {r1:.4f}")
print(f"r0 : {r0:.4f}")
print(f"q(r1) (resonant surface) (pre) : {q_pre_func(r1):.4f}")
print(f"q(r1) (resonant surface) (post) : {q_post_func(r1):.4f}")
print(f"q((r0-r1)/4) (1/4 of the core) (pre) : {q_pre_func(0.25*(r0-r1)):.4f}")
print(f"q((r0-r1)/4) (1/4 of the core) (post) : {q_post_func(0.25*(r0-r1)):.4f}")
print(f"q((r0-r1)/2) (mid core) (pre) : {q_pre_func(0.5*(r0-r1)):.4f}")
print(f"q((r0-r1)/2) (mid core) (post) : {q_post_func(0.5*(r0-r1)):.4f}")
print(f"q(3(r0-r1)/4) (3/4 of the core) (pre) : {q_pre_func(0.75*(r0-r1)):.4f}")
print(f"q(3(r0-r1)/4) (3/4 of the core) (post) : {q_post_func(0.75*(r0-r1)):.4f}")
print(f"q(r0) (boundary radius) (pre) : {q_pre_func(r0):.4f}")
print(f"q(r0) (boundary radius) (post) : {q_post_func(r0):.4f}")