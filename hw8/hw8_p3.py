import numpy as np
from scipy.integrate import quad

np.set_printoptions(precision=8, suppress=True)

print("="*60)
print("Problem 3: Discrete Trigonometric Polynomial S4")
print("f(x) = x^2 * sin(x),  m=16,  interval=[0,1]")
print("="*60)

m = 16
n = 4
N = 2 * m   # 32 data points

# z_i = -pi + (i/m)*pi  =>  x_i = i/(2m)
xi = np.array([i / N for i in range(N)])
zi = np.array([-np.pi + (i / m) * np.pi for i in range(N)])
yi = xi**2 * np.sin(xi)

print(f"\n  2m = {N} data points,  x_i = i/{N},  z_i = -pi + (i/{m})*pi")

# a_k = (1/m) * sum(y_i * cos(k*z_i))
# b_k = (1/m) * sum(y_i * sin(k*z_i)),  k = 1..n-1
a = np.array([(1.0/m) * np.sum(yi * np.cos(k*zi)) for k in range(n+1)])
b = np.zeros(n+1)
for k in range(1, n):
    b[k] = (1.0/m) * np.sum(yi * np.sin(k*zi))

print(f"\n(a) Fourier Coefficients:")
print(f"    a0 = {a[0]:.10f}")
for k in range(1, n):
    print(f"    a{k} = {a[k]:.10f},  b{k} = {b[k]:.10f}")
print(f"    a4 = {a[4]:.10f}")

print(f"\n    S4(z) = (1/2)*a0 + a4*cos(4z)")
print(f"          + a1*cos(z)  + b1*sin(z)")
print(f"          + a2*cos(2z) + b2*sin(2z)")
print(f"          + a3*cos(3z) + b3*sin(3z)")

def S4(z):
    val = 0.5*a[0] + a[n]*np.cos(n*z)
    for k in range(1, n):
        val += a[k]*np.cos(k*z) + b[k]*np.sin(k*z)
    return val

def S4_of_x(xv):
    return S4(np.pi * (2*xv - 1))

# (b) Integral of S4(x) from 0 to 1
# All trig terms vanish => integral = (1/2)*a0
int_S4 = 0.5 * a[0]
int_S4_num, _ = quad(S4_of_x, 0, 1)
print(f"\n(b) Integral of S4(x) from 0 to 1:")
print(f"    = (1/2) * a0 = (1/2) * {a[0]:.10f} = {int_S4:.10f}")
print(f"    scipy.quad verification: {int_S4_num:.10f}")

# (c) Exact integral of x^2*sin(x) from 0 to 1
# Antiderivative: -x^2*cos(x) + 2x*sin(x) + 2*cos(x)
exact = np.cos(1) + 2*np.sin(1) - 2
exact_num, _ = quad(lambda x: x**2 * np.sin(x), 0, 1)
print(f"\n(c) Exact integral of x^2*sin(x) from 0 to 1:")
print(f"    cos(1) + 2*sin(1) - 2 = {exact:.10f}")
print(f"    scipy.quad:             {exact_num:.10f}")
print(f"    |integral_S4 - exact|  = {abs(int_S4 - exact):.2e}")

# (d) Error E(S4) = sum_{i=0}^{2m-1} (y_i - S4(z_i))^2
E_S4 = np.sum((yi - S4(zi))**2)
print(f"\n(d) Error E(S4) = sum(y_i - S4(z_i))^2 = {E_S4:.10f}")

# Formula verification: E = sum(y^2) - m*(0.5*a0^2 + sum_{k=1}^{3}(ak^2+bk^2) + a4^2)
E_formula = np.sum(yi**2) - m*(0.5*a[0]**2 + sum(a[k]**2+b[k]**2 for k in range(1,n)) + a[n]**2)
print(f"    Formula check:          {E_formula:.10f}")
