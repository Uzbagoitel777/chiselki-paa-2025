#Вариант 11.
#Решить систему методом гаусса
#2x1+3x2+3x3+3x4=5
#2x1+2x2+2x3+3x4=5
#2x1+2x2+x3+2x4=4
#2x1+2x2+x3+x4=3
#Предусмотреть верификацию вычисленных корней
from math import isnan

class LinearEquation:
    def __init__(self, coefs:[int|float], value: int|float):
        self.coefs = coefs
        self.value = value
        self.update_nonzeroes()

    def __add__(self, other):
        if len(self) != len(other):
            return AttributeError("Уравнения должны быть одной длины")
        new_eq = LinearEquation(self.coefs, self.value + other.value)
        for i in range(len(self.coefs)):
            new_eq.coefs[i] = self.coefs[i] + other.coefs[i]
        new_eq.update_nonzeroes()
        return new_eq

    def __mul__(self, other: int|float):
        new_eq = LinearEquation(self.coefs.copy(), self.value * other)
        for i in range(len(self)):
            new_eq.coefs[i] = self.coefs[i] * other
        new_eq.update_nonzeroes()
        return new_eq

    def __len__(self):
        return len(self.coefs)

    def __getitem__(self, item):
        return self.coefs[item]

    def __str__(self):
        return self.cast_to_str()

    def update_nonzeroes(self):
        count = 0
        for x in self.coefs:
            try:
                if float(x):
                    count += 1
            except ValueError:
                raise ValueError("Коэффициенты уравнения должны иметь числовой тип данных")
        self.nonzero_coefs_count = count
        return count

    def cast_to_list(self) -> list:
        new_list = []
        for x in self.coefs:
            new_list.append(x)
        return new_list

    def cast_to_str(self) -> str:
        new_str = ''
        for i, coef in enumerate(self.coefs):
            if coef != 0:
                if abs(coef) == 1:
                    new_str += f'x{sup(i + 1, subscript=True)} '
                else:
                    new_str += f'{abs(coef) if i>0 else coef}x{sup(i+1, subscript=True)} '
                if i == len(self)-1:
                    new_str += '= '
                elif self.coefs[i+1] < 0:
                    new_str += '- '
                elif self.coefs[i+1] > 0:
                    new_str += '+ '
            elif i == len(self)-1:
                new_str += '= '
        new_str += str(self.value)
        return new_str

    def calculate(self, roots:[int|float]) -> int|float:
        ans = 0
        for i in range(len(self)):
            ans += self.coefs[i]*roots[i]
        return ans

    def validate_roots(self, roots:[int|float]) -> bool:
        return self.value == self.calculate(roots)

class Matrix:
    def __init__(self, rows: [[int|float]]):
        self.matrix = rows
        if not hasattr(rows[0], 'coefs'):
            self.fill_empty()
        self.row_count = len(rows)
        self.width = len(rows[0])
        self.size = (self.row_count, self.width)

    def fill_empty(self):
        lengths = [len(row) for row in self.matrix]
        max_len = max(lengths)
        for i, row in enumerate(self.matrix):
            if len(row) < max_len:
                self.matrix[i] += [0]*(max_len - len(row))
            self.matrix[i] = [x if not (x is None or isnan(x)) else 0 for x in row]


#System of Linear Equations
class SLE(Matrix):
    def __init__(self, rows: [LinearEquation]):
        super().__init__(rows)
        self.equation_list = rows
        self.values = [row.value for row in rows]
        self.extended = self.update_extended()
        self.source_matrix = rows
        self.source_values = [x for x in self.values]
        self.source_extended = Matrix(self.extended.matrix)
        self.source_equation_list = [x*1 for x in self.equation_list]

    def update_extended(self) -> Matrix:
        self.extended = Matrix([[self.values[i]] + row.coefs for i, row in enumerate(self.matrix)])
        return self.extended

    def pivot_align(self):
        for i in range(self.row_count):
            if self.matrix[i][i] == 0:
                for j in range(i, self.row_count):
                    if self.matrix[j][i] != 0:
                        self.matrix.insert(i, self.matrix.pop(j))

    def gauss_transform(self):
        matrix = self.matrix
        vals = self.values
        for i in range(self.width - 1):
            for j in range(i + 1, self.row_count):
                d = -1*matrix[j][i]/matrix[i][i]
                temp = matrix[i] * d
                matrix[j] += temp
                vals[j] = matrix[j].value
        return self.update_extended()

    def gauss_solve(self) -> [int|float]:
        matrix = self.matrix.copy()
        vals = self.values.copy()
        matrix.reverse()
        vals.reverse()
        last_root = vals[0] / matrix[0][-1]
        roots = [last_root]
        for i in range(1, self.row_count):
            root = vals[i]
            for j in range(1, i + 1):
                root -= matrix[i][-j] * roots[j - 1]
            roots.append(root/matrix[i][-(i+1)])

        return roots[::-1]


    def verify_solution(self) -> bool:
        self.gauss_transform()
        print('СЛАУ после преобразования:')
        for row in self.equation_list:
            print(row)
        solution = self.gauss_solve()
        print('Корни, найденные методом Гаусса:')
        for i, root in enumerate(solution):
            print(f'x{sup(i, subscript=True)} = {root}')
        valid = True
        for i, equation in enumerate(self.source_equation_list):
            if not equation.validate_roots(solution):
                valid = False
            print(f'''
Проверка уравнения #{i}: {str(equation)}
{equation.validate_roots(solution)}
Ожидалось: {equation.value}. Получено: {equation.calculate(solution)}''')
        print(f"Итого, правильно ли решение? {valid}")
        return valid

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

#2x1+3x2+3x3+3x4=5
#2x1+2x2+2x3+3x4=5
#2x1+2x2+x3+2x4=4
#2x1+2x2+x3+x4=3
def main():
    eqs = [
        LinearEquation([2, 3, 3, 3], 5),
        LinearEquation([2, 2, 2, 3], 5),
        LinearEquation([2, 2, 1, 2], 4),
        LinearEquation([2, 2, 1, 1], 3)
    ]
    print('СЛАУ:')
    for eq in eqs:
        print(eq)
    sle = SLE(eqs)
    sle.verify_solution()

if __name__ == '__main__':
    main()






