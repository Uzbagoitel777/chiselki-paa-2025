import matplotlib.pyplot as plt
import numpy as np

def linear_approximation(x, y):
    x = np.array(x)
    y = np.array(y)
    n = len(x)

    A = np.array([
        [np.sum(x * x), np.sum(x)],
        [np.sum(x), n]
    ])

    B = np.array([
        np.sum(x * y),
        np.sum(y)
    ])

    a, b = map(float, np.linalg.solve(A, B))
    return a, b


def quadratic_approximation(x, y):
    x = np.array(x)
    y = np.array(y)
    n = len(x)

    A = np.array([
        [np.sum(x**4), np.sum(x**3), np.sum(x**2)],
        [np.sum(x**3), np.sum(x**2), np.sum(x)],
        [np.sum(x**2), np.sum(x), n]
    ])

    B = np.array([
        np.sum(x**2 * y),
        np.sum(x * y),
        np.sum(y)
    ])

    a, b, c = map(float, np.linalg.solve(A, B))
    return a, b, c


def linear_value(x, a, b):
    return a * x + b


def quadratic_value(x, a, b, c):
    return a * x * x + b * x + c


def plot_approximation(x_data, y_data, coeffs, kind):
    x_data = np.array(x_data)
    y_data = np.array(y_data)

    x_min = x_data.min()
    x_max = x_data.max()
    x_plot = np.linspace(x_min, x_max, 500)

    if kind == "linear":
        a, b = coeffs
        y_plot = linear_value(x_plot, a, b)
    elif kind == "quadratic":
        a, b, c = coeffs
        y_plot = quadratic_value(x_plot, a, b, c)
    else:
        return

    plt.scatter(x_data, y_data, marker='o', alpha=0.95, color='orange')
    plt.plot(x_plot, y_plot,  linestyle='-', linewidth=2, color='red')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()


def main():
    # x1, y1 = ([0, 1, 2, 4], [0.2, 0.9, 2.1, 3.7])
    # x2, y2 = ([-2, -1, 0, 1, 2], [6, 2, -1, -2, -1])
    # ex_linear = linear_approximation(x1, y1)
    # ex_quadratic = quadratic_approximation(x2, y2)
    # print(ex_linear)
    # print(ex_quadratic)
    # plot_approximation(x1, y1, ex_linear, 'linear')
    # plot_approximation(x2, y2, ex_quadratic, 'quadratic')

    x = [-1, 2, 3, 4]
    y = [1, 7, 9, 11]
    linear = linear_approximation(x, y)
    quadratic = quadratic_approximation(x, y)
    a, b = linear
    print(f'Линейная аппроксимация: {a}x {"-" if b < 0 else "+"} {abs(b)}')
    a, b, c = quadratic
    print(f'Полиномиальная аппроксимация 2-й степени: {a}x² {"-" if b < 0 else "+"} {abs(b)}x {"-" if c < 0 else "+"} {abs(c)}')
    plot_approximation(x, y, linear, 'linear')
    plot_approximation(x, y, quadratic, 'quadratic')


if __name__ == '__main__':
    main()