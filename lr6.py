import numpy as np
import matplotlib.pyplot as plt


def linear_approximation(x, y):
    A1 = np.vstack([x, np.ones(len(x))]).T
    a1, b1 = np.linalg.lstsq(A1, y, rcond=None)[0]
    return a1, b1


def polynomial_approximation(x, y):
    A2 = np.vstack([x ** 2, x, np.ones(len(x))]).T
    a2, b2, c2 = np.linalg.lstsq(A2, y, rcond=None)[0]
    return a2, b2, c2


def evaluate_linear(x, a, b):
    return a * x + b


def evaluate_polynomial(x, a, b, c):
    return a * x ** 2 + b * x + c


def calculate_errors(y_actual, y_predicted):
    mse = np.mean((y_actual - y_predicted) ** 2)
    rmse = np.sqrt(mse)
    return mse, rmse


if __name__ == '__main__':
    x = np.array([-1, 2, 3, 4])
    y = np.array([1, 7, 9, 11])

    print("Исходные данные:")
    print("X:", x)
    print("Y:", y)
    print()

    a1, b1 = linear_approximation(x, y)

    print("ЛИНЕЙНАЯ АППРОКСИМАЦИЯ")
    print(f"Формула: y = a·x + b")
    print(f"\nКоэффициенты:")
    print(f"  a = {a1:.4f}")
    print(f"  b = {b1:.4f}")
    print(f"\nИтоговая зависимость: y = {a1:.4f}·x + {b1:.4f}")
    print()

    a2, b2, c2 = polynomial_approximation(x, y)

    print("ПОЛИНОМИАЛЬНАЯ АППРОКСИМАЦИЯ (2-я степень)")
    print(f"Формула: y = a·x² + b·x + c")
    print(f"\nКоэффициенты:")
    print(f"  a = {a2:.4e}")
    print(f"  b = {b2:.4f}")
    print(f"  c = {c2:.4f}")
    print(f"\nИтоговая зависимость: y = {a2:.4e}·x² + {b2:.4f}·x + {c2:.4f}")
