import math


def simpson_weights(count):
    if count % 2 != 0:
        raise ValueError("n and m must be even.")
    weights = [0] * (count + 1)
    weights[0] = 1
    weights[-1] = 1
    for i in range(1, count):
        weights[i] = 4 if i % 2 else 2
    return weights


def double_simpson(func, left, right, low, high, n, m):
    hx = (right - left) / (2 * n)
    wx = simpson_weights(2 * n)
    total = 0.0

    for i in range(2 * n + 1):
        x = left + i * hx
        yl = low(x)
        yu = high(x)
        hy = (yu - yl) / (2 * m)
        wy = simpson_weights(2 * m)
        inner_sum = 0.0
        for j in range(2 * m + 1):
            y = yl + j * hy
            inner_sum += wy[j] * func(x, y)
        total += wx[i] * (hy / 3.0) * inner_sum

    return total * hx / 3.0


def gauss_table(n):
    if n == 3:
        return (
            [-0.775, 0.0, 0.775],
            [0.556, 0.889, 0.556],
        )
    raise ValueError("Only n=3 is supported.")


def double_gauss(func, left, right, low, high, n, m):
    x_nodes, x_weights = gauss_table(n)
    y_nodes, y_weights = gauss_table(m)
    total = 0.0

    for x_hat, wx in zip(x_nodes, x_weights):
        x = 0.5 * (right - left) * x_hat + 0.5 * (left + right)
        jac_x = 0.5 * (right - left)
        yl = low(x)
        yu = high(x)
        jac_y = 0.5 * (yu - yl)

        for y_hat, wy in zip(y_nodes, y_weights):
            y = jac_y * y_hat + 0.5 * (yl + yu)
            total += wx * wy * func(x, y) * jac_x * jac_y

    return total


def main():
    f = lambda x, y: 2.0 * y * math.sin(x) + math.cos(x) ** 2
    left, right = 0.0, math.pi / 4.0
    low = lambda x: math.sin(x)
    high = lambda x: math.cos(x)

    simpson_nm = double_simpson(f, left, right, low, high, n=4, m=4)
    gauss_nm = double_gauss(f, left, right, low, high, n=3, m=3)
    reference = double_simpson(f, left, right, low, high, n=200, m=200)

    print("HW4 Q3: double integral with variable y-bounds")
    print(f"Simpson (n=4,m=4) : {simpson_nm:.12f} (error={abs(simpson_nm - reference):.3e})")
    print(f"Gaussian (n=3,m=3): {gauss_nm:.12f} (error={abs(gauss_nm - reference):.3e})")
    print(f"Reference value    : {reference:.12f}")


if __name__ == "__main__":
    main()
