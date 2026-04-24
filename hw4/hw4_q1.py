import math


def trap_rule(func, left, right, step):
    count = int(round((right - left) / step))
    grid = [left + i * step for i in range(count + 1)]
    return step * (0.5 * func(grid[0]) + sum(func(x) for x in grid[1:-1]) + 0.5 * func(grid[-1]))


def midpoint_rule(func, left, right, step):
    two_n = int(round((right - left) / step))
    if two_n % 2 != 0:
        raise ValueError("(b-a)/h must be even.")
    half_n = two_n // 2
    return 2.0 * step * sum(func(left + (2 * i - 1) * step) for i in range(1, half_n + 1))


def simpson_rule(func, left, right, step):
    count = int(round((right - left) / step))
    if count % 2 != 0:
        raise ValueError("Need even subintervals for Simpson.")
    grid = [left + i * step for i in range(count + 1)]
    odd = sum(func(grid[i]) for i in range(1, count, 2))
    even = sum(func(grid[i]) for i in range(2, count, 2))
    return step / 3.0 * (func(grid[0]) + 4.0 * odd + 2.0 * even + func(grid[-1]))


def exact_value(left, right):
    def antiderivative(x):
        return math.exp(x) * (math.sin(4.0 * x) - 4.0 * math.cos(4.0 * x)) / 17.0

    return antiderivative(right) - antiderivative(left)


def main():
    left, right, step = 1.0, 2.0, 0.1
    func = lambda x: math.exp(x) * math.sin(4.0 * x)

    trap = trap_rule(func, left, right, step)
    simp = simpson_rule(func, left, right, step)
    mid = midpoint_rule(func, left, right, step)
    exact = exact_value(left, right)

    print("HW4 Q1: ∫[1,2] e^x sin(4x) dx, h=0.1")
    print(f"Composite Trapezoidal : {trap:.12f} (error={abs(trap - exact):.3e})")
    print(f"Composite Simpson     : {simp:.12f} (error={abs(simp - exact):.3e})")
    print(f"Composite Midpoint    : {mid:.12f} (error={abs(mid - exact):.3e})")
    print(f"Exact value           : {exact:.12f}")


if __name__ == "__main__":
    main()
