import numpy as np

np.set_printoptions(precision=8, suppress=True)

x = np.array([4.0, 4.2, 4.5, 4.7, 5.1, 5.5, 5.9, 6.3])
y = np.array([102.6, 113.2, 130.1, 142.1, 167.5, 195.1, 224.9, 256.8])
n = len(x)

print("="*60)
print("Problem 1: Least Squares Approximations from Data")
print("="*60)
print(f"Data: x = {x}")
print(f"      y = {y}")

# ---- 1a: Degree-2 polynomial P(x) = c0 + c1*x + c2*x^2 ----
print("\n--- 1a: Polynomial Degree 2: P(x) = c0 + c1*x + c2*x^2 ---")
A1a = np.column_stack([np.ones(n), x, x**2])
c1a = np.linalg.solve(A1a.T @ A1a, A1a.T @ y)
print(f"  c0 = {c1a[0]:.6f}")
print(f"  c1 = {c1a[1]:.6f}")
print(f"  c2 = {c1a[2]:.6f}")
print(f"  P(x) = {c1a[0]:.6f} + {c1a[1]:.6f}*x + {c1a[2]:.6f}*x^2")
y1a = A1a @ c1a
E1a = np.sum((y - y1a)**2)
print(f"  Error E = sum(y_i - P(x_i))^2 = {E1a:.6f}")

print("\n  Verification (x, y_actual, y_fitted, residual):")
for i in range(n):
    print(f"    x={x[i]:.1f}: y={y[i]:.1f},  P(x)={y1a[i]:.4f},  res={y[i]-y1a[i]:.4f}")

# ---- 1b: y = b*exp(a*x) ----
print("\n--- 1b: Exponential: y = b * exp(a*x) ---")
lnY = np.log(y)
A1b = np.column_stack([np.ones(n), x])
c1b = np.linalg.solve(A1b.T @ A1b, A1b.T @ lnY)
b1b = np.exp(c1b[0])
a1b = c1b[1]
print(f"  b = {b1b:.6f}")
print(f"  a = {a1b:.6f}")
print(f"  y = {b1b:.6f} * exp({a1b:.6f} * x)")
y1b = b1b * np.exp(a1b * x)
E1b = np.sum((y - y1b)**2)
print(f"  Error E = sum(y_i - b*exp(a*x_i))^2 = {E1b:.6f}")

print("\n  Verification (x, y_actual, y_fitted, residual):")
for i in range(n):
    print(f"    x={x[i]:.1f}: y={y[i]:.1f},  fit={y1b[i]:.4f},  res={y[i]-y1b[i]:.4f}")

# ---- 1c: y = b*x^n ----
print("\n--- 1c: Power law: y = b * x^n ---")
A1c = np.column_stack([np.ones(n), np.log(x)])
c1c = np.linalg.solve(A1c.T @ A1c, A1c.T @ lnY)
b1c = np.exp(c1c[0])
n1c = c1c[1]
print(f"  b = {b1c:.6f}")
print(f"  n = {n1c:.6f}")
print(f"  y = {b1c:.6f} * x^{n1c:.6f}")
y1c = b1c * x**n1c
E1c = np.sum((y - y1c)**2)
print(f"  Error E = sum(y_i - b*x_i^n)^2 = {E1c:.6f}")

print("\n  Verification (x, y_actual, y_fitted, residual):")
for i in range(n):
    print(f"    x={x[i]:.1f}: y={y[i]:.1f},  fit={y1c[i]:.4f},  res={y[i]-y1c[i]:.4f}")
