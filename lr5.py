import numpy as np
import matplotlib.pyplot as plt


def calculate_cubic_spline(x, y):
    n = len(x) - 1

    h = [x[i + 1] - x[i] for i in range(n)]
    alpha = [0.0] * (n)
    for i in range(1, n):
        alpha[i] = (3 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1])

    l = [1.0] + [0.0] * n
    mu = [0.0] * (n + 1)
    z = [0.0] * (n + 1)
    for i in range(1, n):
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    l[n] = 1.0
    z[n] = 0.0
    c = [0.0] * (n + 1)
    b = [0.0] * n
    d = [0.0] * n

    for j in range(n - 1, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])

    coeffs = []
    for i in range(n):
        coeffs.append([y[i], b[i], c[i], d[i]])

    return coeffs


def spline_value(x_val, x, coeffs):
    n = len(x) - 1
    if x_val < x[0] or x_val > x[-1]:
        print("Внимание: точка вне диапазона исходных данных!")
    i = 0
    while i < n and x_val > x[i + 1]:
        i += 1
    if i == n:
        i = n - 1
    dx = x_val - x[i]
    a, b_coef, c_coef, d_coef = coeffs[i]
    return a + b_coef * dx + c_coef * dx ** 2 + d_coef * dx ** 3


def print_calculation_process(x_input, x, coeffs, y_output):
    n = len(x) - 1
    interval = 0
    while interval < n and x_input > x[interval + 1]:
        interval += 1
    if x_input == x[-1]:
        interval = n - 1

    print(f"1. X = {x_input} находится в интервале [{x[interval]}, {x[interval + 1]}]")
    print(f"2. Коэффициенты сплайна на этом интервале:")
    print(f"   a (y[{interval}]) = {coeffs[interval][0]}")
    print(f"   b = {coeffs[interval][1]:.4f}")
    print(f"   c = {coeffs[interval][2]:.4f}")
    print(f"   d = {coeffs[interval][3]:.4f}")
    dx = x_input - x[interval]
    print(f"3. Вычисляем dx = X - x[{interval}] = {dx:.4f}")
    print(f"4. Y = a + b*dx + c*dx² + d*dx³")
    print(f"   = {coeffs[interval][0]:.4f} + {coeffs[interval][1]:.4f}*{dx:.4f} + " +
          f"{coeffs[interval][2]:.4f}*{dx ** 2:.4f} + {coeffs[interval][3]:.4f}*{dx ** 3:.4f}")
    print(f"   = {y_output:.4f}")


def plot_spline(x, y, coeffs, x_input, y_output):
    x_plot = np.linspace(x[0], x[-1], 200)
    y_plot = [spline_value(xi, x, coeffs) for xi in x_plot]

    plt.figure(figsize=(10, 6))
    plt.plot(x_plot, y_plot, 'b-', label='Кубический сплайн', linewidth=2)
    plt.plot(x, y, 'ro', label='Исходные точки', markersize=8)
    plt.scatter([x_input], [y_output], color='green', s=150, zorder=5,
                label=f'Точка запроса ({x_input:.2f}, {y_output:.2f})')
    plt.xlabel('X', fontsize=12)
    plt.ylabel('Y', fontsize=12)
    plt.title('Кубический сплайн для функции №21', fontsize=14)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    x = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    y = [2.8, 2.0, 1.8, 1.6, 2.2, 3.0]

    coeffs = calculate_cubic_spline(x, y)

    print("Исходные точки:")
    for i in range(len(x)):
        print(f"x={x[i]}, y={y[i]}")
    print()

    x_input = float(input(f"Введите X в интервале {x[0]} - {x[-1]}: "))
    y_output = spline_value(x_input, x, coeffs)
    print(f"\nПри X = {x_input}")
    print(f"Y = {y_output:.4f}")

    print("\nПроцесс вычисления:")
    print_calculation_process(x_input, x, coeffs, y_output)

    plot_spline(x, y, coeffs, x_input, y_output)