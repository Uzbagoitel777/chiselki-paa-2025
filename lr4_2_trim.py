import matplotlib.pyplot as plt

def lagrange_interpolation(x_vals, y_vals):
    n = len(x_vals)
    result = [0.0] * n
    for i in range(n):
        basis = [y_vals[i]]
        for j in range(n):
            if i != j:
                denominator = x_vals[i] - x_vals[j]

                new_basis = [0.0] * (len(basis) + 1)
                for k in range(len(basis)):
                    new_basis[k] += basis[k] * (-x_vals[j])  # constant term
                    new_basis[k + 1] += basis[k]  # x term

                basis = [coef / denominator for coef in new_basis]
        for k in range(len(basis)):
            result[k] += basis[k]

    return result


def newton_interpolation(x_vals, y_vals):
    n = len(x_vals)
    divided_diff = [[0.0] * n for _ in range(n)]
    for i in range(n):
        divided_diff[i][0] = y_vals[i]
    for j in range(1, n):
        for i in range(n - j):
            divided_diff[i][j] = (divided_diff[i + 1][j - 1] - divided_diff[i][j - 1]) / (x_vals[i + j] - x_vals[i])

    newton_coeffs = [divided_diff[0][i] for i in range(n)]
    result = [0.0] * n
    result[0] = newton_coeffs[0]
    current_poly = [1.0]
    for i in range(1, n):
        new_poly = [0.0] * (len(current_poly) + 1)
        for j in range(len(current_poly)):
            new_poly[j] += current_poly[j] * (-x_vals[i - 1])  # constant term
            new_poly[j + 1] += current_poly[j]  # x term

        current_poly = new_poly
        for j in range(len(current_poly)):
            result[j] += newton_coeffs[i] * current_poly[j]

    return result


def evaluate_polynomial(coeffs, x):
    result = 0.0
    for i, coef in enumerate(coeffs):
        result += coef * (x ** i)
    return result


def plot_interpolations(x_vals, y_vals, num_points=200):
    lagrange_coeffs = lagrange_interpolation(x_vals, y_vals)
    newton_coeffs = newton_interpolation(x_vals, y_vals)

    x_min, x_max = min(x_vals), max(x_vals)
    x_range = x_max - x_min
    x_plot = [x_min - 0.1 * x_range + i * (1.2 * x_range) / (num_points - 1)
              for i in range(num_points)]

    y_lagrange = [evaluate_polynomial(lagrange_coeffs, x) for x in x_plot]
    y_newton = [evaluate_polynomial(newton_coeffs, x) for x in x_plot]

    plt.figure(figsize=(10, 6))
    plt.plot(x_plot, y_lagrange, 'b-', label='Лагранж', linewidth=2)
    plt.plot(x_plot, y_newton, 'r--', label='Ньютон', linewidth=2, alpha=0.7)
    plt.scatter(x_vals, y_vals, color='black', s=100, zorder=5, label='')

    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.title('График интерполированных многочленов', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    x = [-2, 0, 1]
    y = [4, 1, 3]
    lagr_ex = lagrange_interpolation(x, y)
    newton_ex = newton_interpolation(x, y)
    print("Решение варианта Лагранж:", lagr_ex)
    print("Решение варианта Ньютон:", newton_ex)
    print(f"Лагранж при х = 10: {evaluate_polynomial(lagr_ex, 10)}")
    print(f"Ньютон при х = 10: {evaluate_polynomial(newton_ex, 10)}")
    plot_interpolations(x, y)

