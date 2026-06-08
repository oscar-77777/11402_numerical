import numpy as np
from scipy.integrate import quad

np.set_printoptions(precision=8, suppress=True)

print("="*60)
print("Problem 2: Continuous Least Squares Polynomial on [-1,1]")
print("f(x) = (1/2)*cos(x) + (1/4)*sin(2x)")
print("="*60)

def f(x):
    return 0.5 * np.cos(x) + 0.25 * np.sin(2*x)

phi = [lambda x: np.ones_like(np.atleast_1d(x)),
       lambda x: np.atleast_1d(x).astype(float),
       lambda x: np.atleast_1d(x).astype(float)**2]

# Gram matrix G[i,j] = integral_{-1}^{1} phi_i(x)*phi_j(x) dx
G = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        G[i, j], _ = quad(lambda x, i=i, j=j: phi[i](x) * phi[j](x), -1, 1)

# RHS b[i] = integral_{-1}^{1} f(x)*phi_i(x) dx
rhs = np.zeros(3)
for i in range(3):
    rhs[i], _ = quad(lambda x, i=i: f(x) * phi[i](x), -1, 1)

print(f"\n  Gram matrix G:")
for row in G:
    print(f"    {row}")
print(f"\n  RHS b = {rhs}")

a = np.linalg.solve(G, rhs)
print(f"\n  a0 = {a[0]:.8f}")
print(f"  a1 = {a[1]:.8f}")
print(f"  a2 = {a[2]:.8f}")
print(f"  P2(x) = {a[0]:.6f} + {a[1]:.6f}*x + {a[2]:.6f}*x^2")

def P2(x):
    return a[0] + a[1]*x + a[2]*x**2

E, _ = quad(lambda x: (f(x) - P2(x))**2, -1, 1)
print(f"  Error E = integral_(-1)^(1) (f - P2)^2 dx = {E:.10f}")

print("\n  Verification at sample points:")
xs = np.linspace(-1, 1, 9)
for xi in xs:
    print(f"    x={xi:+.3f}: f(x)={f(xi):.6f},  P2(x)={P2(xi):.6f},  diff={f(xi)-P2(xi):.6f}")
