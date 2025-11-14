#Решить систему уравнений методом прогонки
#10x1-4x2=8
#x1+2x2-0.2x3=5.5
#x2-7x3+x4=2
#-2x3+5x4=-1
from lec3 import SLE, LinearEquation, sup

def diagonal_solve(datamatrix: [[float]]) -> [float]:
    size = len(datamatrix)
    u_arr, v_arr = [], []
    a, b, c, d = datamatrix[0]
    u_arr.append(-1*c/b)
    v_arr.append(d/b)
    for i in range(1, size):
        a, b, c, d = datamatrix[i]
        u_prev, v_prev = u_arr[i-1], v_arr[i-1]
        denom = a*u_prev + b
        u_arr.append(-1*c / denom)
        v_arr.append((d - a*v_prev) / denom)
    x_arr = [0]*size
    x_arr[-1] = v_arr[-1]
    for i in range(size-2, -1, -1):
        x_arr[i] = (u_arr[i]*x_arr[i+1] + v_arr[i])

    return x_arr


def main():
    # [a, b, c, d]

    # Пример из лекции
    mtx = [
        [0, 10, 1, 5],
        [-2, 9, 1, -1],
        [0.1, 4, -1, -5],
        [-1, 8, 0, 40]
    ]
    # Вариант 11
    mtx_2 = [
        [0, 10, 4, 8],
        [1, 2, -0.2, 5.5],
        [1, -7, 1, 2],
        [-2, 5, 0, -1]
    ]
    print(diagonal_solve(mtx))
    print(diagonal_solve(mtx_2))


if __name__ == '__main__':
    main()
