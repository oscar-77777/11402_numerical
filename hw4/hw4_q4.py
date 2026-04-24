import math


def composite_simpson(f, a, b, n):
    if n % 2 != 0:
        raise ValueError("Simpson's rule requires even n.")
    h = (b - a) / n
    s = f(a) + f(b)
    s += 4.0 * sum(f(a + i * h) for i in range(1, n, 2))
    s += 2.0 * sum(f(a + i * h) for i in range(2, n, 2))
    return s * h / 3.0


def fa(t):
    return -(t ** (-7.0 / 4.0)) * math.sin(1.0 / t)


def fb(t):
    return -(t ** 2) * math.sin(1.0 / t)


def main():
    n = 4
    i1 = composite_simpson(fa, 9.0, 1.0, 2 * n)
    i2 = composite_simpson(fb, 1.0, 1e-4, 2 * n)

    print("HW4 Q4")
    print(f"a: {i1:.6f}")
    print(f"b: {i2:.6f}")


if __name__ == "__main__":
    main()
