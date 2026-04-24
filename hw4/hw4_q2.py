import math


def gauss_quad(func, left, right, n):
    coeff = {
        3: (
            [-0.775, 0.0, 0.775],
            [0.556, 0.889, 0.556],
        ),
        4: (
            [
                -0.861,
                -0.340,
                0.340,
                0.861,
            ],
            [
                0.348,
                0.652,
                0.652,
                0.348,
            ],
        ),
    }
    if n not in coeff:
        raise ValueError("n must be 3 or 4.")

    nodes, weights = coeff[n]
    mid = 0.5 * (left + right)
    half = 0.5 * (right - left)
    return half * sum(w * func(mid + half * x) for x, w in zip(nodes, weights))


def exact_value(left, right):
    def antiderivative(x):
        return (x ** 3 / 3.0) * math.log(x) - x ** 3 / 9.0

    return antiderivative(right) - antiderivative(left)


def main():
    left, right = 1.0, 1.5
    func = lambda x: (x ** 2) * math.log(x)

    g3 = gauss_quad(func, left, right, 3)
    g4 = gauss_quad(func, left, right, 4)
    exact = exact_value(left, right)

    print("HW4 Q2: ∫[1,1.5] x^2 ln(x) dx")
    print(f"Gaussian Quadrature n=3: {g3:.12f} (error={abs(g3 - exact):.3e})")
    print(f"Gaussian Quadrature n=4: {g4:.12f} (error={abs(g4 - exact):.3e})")
    print(f"Exact value            : {exact:.12f}")


if __name__ == "__main__":
    main()
