import matplotlib.pyplot as plt
import numpy as np

class Polynomial:
    def __init__(self, max_pow: int, coefs: [float | int]):
        self.max_pow = max_pow
        self.coefs = coefs
        self.len = len(coefs)
        if self.len < max_pow:
            self.coefs.append([0]*(max_pow - self.len))

    def __add__(self, other):
        if type(other) is float or type(other) is int:
            new_coefs = self.coefs.copy()
            new_coefs[0] += other
            return Polynomial(self.max_pow, new_coefs)
        new_max_pow = max(self.max_pow, other.max_pow)
        min_pow = min(self.max_pow, other.max_pow)
        new_poly = Polynomial(new_max_pow, [0]*(new_max_pow+1))
        for i in range(min_pow+1):
            new_poly.coefs[i] = self.coefs[i] + other.coefs[i]
        return new_poly

    def __sub__(self, other):
        inverted = Polynomial(other.max_pow, [-x for x in other.coefs])
        return self + inverted

    # Только умножение многочлена на константу
    def __mul__(self, other: float | int):
        return Polynomial(self.max_pow, [x*other for x in self.coefs])

    def compute(self, x=1) -> float:
        out = 0
        for i, coef in enumerate(self.coefs):
            out += coef * x**i
        return out

    def to_string(self, show_zero = False) -> str:
        f_coefs = self.coefs.copy()
        terms = ['']*(self.max_pow+1)
        terms[0] = str(f_coefs.pop(0))
        if terms[0] == '0' and not show_zero:
            terms[0] = ''
        for i, coef in enumerate(f_coefs, 1):
            if coef != 0 and not show_zero:
                if coef != 1:
                    terms[i] = str(coef)
                terms[i] += 'x'
                if i > 1:
                    terms[i] += sup(i)
                    if coef > 0:
                        terms[i] += '+'
                    else:
                        terms[i] += '-'

        out = ''
        for term in list(reversed(terms)):
            out += term
        return out



def poly_derivative(poly: Polynomial, depth=1, no_warning = False) -> Polynomial:
    new_max_pow = poly.max_pow - 1
    if new_max_pow < 0:
        if not no_warning:
            print('''Warning: the polynomial is a constant or does not have a derivative. Returning 0
    Set show_zero=True in Polynomial.to_string() to display it
    To silence this warning set no_warning=True in poly_derivative()''')
        return Polynomial(0, [0])
    new_coefs = [0]*(new_max_pow+1)
    for i, coef in enumerate(poly.coefs[1:], 1):
        new_coefs[i - 1] = i * coef
    new_poly = Polynomial(new_max_pow, new_coefs)
    if depth == 1:
        return new_poly
    else:
        return poly_derivative(new_poly, depth = depth-1)

superscripts = {
    '1': '¹', '2': '²', '3': '³', '4': '⁴',
    '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸',
    '9': '⁹', 'n': 'ⁿ', 'i': 'ⁱ', '0': '⁰'
}
subscripts = {
    '1': '₁', '2': '₂', '3': '₃', '4': '₄',
    '5': '₅', '6': '₆', '7': '₇', '8': '₈',
    '9': '₉', 'n': 'ₙ', 'k': 'ₖ', '0': '₀'
}


def sup(x, subscript=False):
    x = str(x)
    out = ''
    for i in range(len(x)):
        if subscript:
            out = out + subscripts[x[i]]
        else:
            out = out + superscripts[x[i]]

    return out


def approx_newton(poly: Polynomial, interval: (float, float), epsilon=0.01, to_zero=False):
    print(f'''Начинаем искать приближённый корень уравнения {poly.to_string()}=0 на интервале {interval} методом Ньютона 
    Точность \u03b5: {epsilon}''')
    delta_x = 1000
    first_derivative = poly_derivative(poly, no_warning=True)
    second_derivative = poly_derivative(poly, depth=2, no_warning=True)
    #Если F(a)F''(a) > 0, то x0=a, иначе b, где a и b - начало и конец отрезка
    x0 = interval[0] if poly.compute(interval[0]) * second_derivative.compute(interval[0]) > 0 else interval[1]
    print(f"Т.к. F({x0})F''({x0}) > 0, x{subscripts['0']} = {x0}")
    iterations = [(poly, x0)]
    x = x0
    while delta_x > epsilon:
        i = len(iterations)
        # y = (x-xk)F'(xk) + F(xk)
        iteration = (Polynomial(1, [x*-1, 1])*first_derivative.compute()
                     + poly.compute(x))
        print(f'Приближение #{i}: {iteration.to_string()}')
        prev_x = iterations[i-1][1]
        x = prev_x - (poly.compute(prev_x) / first_derivative.compute(prev_x))
        print(f'x{sup(i, True)} = {x}')
        iterations.append((iteration, x))
        if to_zero:
            delta_x = abs(poly.compute(x))
            print(f'|F({x})| = {delta_x}')
        else:
            delta_x = abs(x - prev_x)
            print(f'Корень находится в интервале [{min(x, prev_x)}, {max(x, prev_x)}]')
            print(f'|x{sup(i, True)} - x{sup(i-1, True)}| = {delta_x}')
    print(f'''Т.к. {f"F(x{subscripts['k']})" if to_zero else "\u0394x"} = {delta_x} < \u03b5 = {epsilon},
можем остановиться на {i}-м приближении.
Ответ: 
    Приближённым корнем уравнения {poly.to_string()}=0 является x={round(x, 2)}''')
    return {'result': x, 'iterations': iterations}


def approx_binary(poly: Polynomial, interval: (float, float), epsilon=0.01):
    print(f'''Начинаем искать приближённый корень уравнения {poly.to_string()}=0 на интервале {interval} методом деления отрезка пополам 
        Точность \u03b5: {epsilon}''')
    delta_x = float('inf')
    x0 = sum(interval)/2
    iterations = [(interval, x0)]
    x = x0
    while delta_x > epsilon:
        i = len(iterations)
        prev_interval = iterations[i-1][0]
        for bound in prev_interval:
            if poly.compute(x) * poly.compute(bound) < 0:
                new_bound = bound
        new_interval = (x, new_bound)
        x = sum(new_interval)/2
        iterations.append((new_interval, x))
        delta_x = abs(max(new_interval) - min(new_interval))
        print(f'Корень находится в интервале {new_interval}')
    print(f'''Т.к. длина нового отрезка = {delta_x} < \u03b5 = {epsilon},
    можем остановиться на {i}-м приближении.
    Ответ: 
        Приближённым корнем уравнения {poly.to_string()}=0 является x={round(x, 2)}''')
    return {'result': x, 'iterations': iterations}

def graph_approx(poly: Polynomial, iterations):
    x = np.linspace(-1, 2, 100)

    fig = plt.figure(figsize=(14, 8))

    y = poly.compute(x)
    plt.plot(x, y, label=poly.to_string())
    y1 = x*0
    plt.plot(x, y1, label='0')
    for i in range(1, len(iterations)):
        iter = iterations[i][0]
        y = iter.compute(x)
        if i < len(iterations)-1:
            plt.plot(x, y, label=f'Приближение №{i}: {iter.to_string()}')
        else:
            plt.plot(x, y, 'r', label=f'Приближение №{i}: {iter.to_string()}')

    plt.legend()
    plt.grid(True, linestyle=':')
    plt.xlim([0.93, 0.95])
    plt.ylim([-0.02, 0.05])
    plt.title('График уравнения и промежуточных касательных')
    plt.xlabel('ось x')
    plt.ylabel('ось y')
    plt.show()

def test():
    f1 = Polynomial(3, [-1, 1, 0, 1])
    x = 2
    print(f'{f1.to_string()} at x = {x} equals {f1.compute()}')
    print(f'First derivative of {f1.to_string()} is {poly_derivative(f1).to_string()}')
    print(f'Second derivative of {f1.to_string()} is {poly_derivative(f1, 2).to_string()}')
    print(f'Third derivative of {f1.to_string()} is {poly_derivative(f1, 3).to_string()}')
    # print(f'Fourth derivative of {f1.to_string()} is {poly_derivative(f1, 4).to_string(show_zero=True)}')

    approx_newton(f1, (0, 1))


def main():
    f = Polynomial(3, [-1.2, 0.4, 0, 1])
    iters = approx_newton(f, (0, 1))['iterations']
    approx_binary(f, (0,1))
    graph_approx(f, iters)

if __name__ == '__main__':
    main()

