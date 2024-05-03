import math

from tabulate import tabulate
import os.path
from sympy import diff, symbols, cos, sin, sqrt, lambdify, sympify, tan, exp, pi, solve, ln
from check import typeofread, readabe, is_number
from Matrix import M
import matplotlib.pyplot as plt
import numpy as np

x, y = symbols('x y')


def printf(stri, fi, file):
    if fi:
        file.write(str(stri) + "\n")
    else:
        print(stri)


def fun(fr, x1, y1, a, fi2, f2, R2, j):
    printf("Коэффициенты аппроксимирующей функции " + str(fr) + ":", fi2, f2)
    for i in range(len(a)):
        if fi2:
            f2.write(str(a[i]) + " ")
        else:
            print(a[i], end=" ")
    printf("", fi2, f2)
    frc = lambdify(x, fr)
    rez = []
    si = 0
    phi = [0] * n
    for i in range(n):
        phi[i] = frc(x1[i])
        si += (frc(x1[i]) - y1[i]) ** 2
        rez.append([x1[i], y1[i], frc(x1[i]), frc(x1[i]) - y1[i]])
    printf("Мера отклонения:" + str(si), fi2, f2)
    si = (si / n) ** (1 / 2)
    ch, zn = 0, 0
    for i in range(n):
        ch += (y1[i] - phi[i]) ** 2
        zn += (y1[i] - sum(phi) / n) ** 2
    R2[j] = 1 - ch / zn

    printf("Среднеквадратичное отклонение:" + str(si), fi2, f2)
    printf("Коэффициент детерминации: " + str(R2[j]), fi2, f2)
    if R2[j] >= 0.95:
        printf("модель хорошо описывает явление", fi2, f2)
    elif R2[j] >= 0.75:
        printf("модель в целом адекватно описывает явление", fi2, f2)
    elif R2[j] >= 0.5:
        printf("модель слабо описывает явление", fi2, f2)
    else:
        printf("точность аппроксимации недостаточна и модель требует изменения", fi2, f2)
    printf(tabulate(rez, ["x", "y", "phi(x)", "e"]), fi2, f2)
    printf("", fi2, f2)
    return R2


u = typeofread(False, "", "Выберете откуда вводите информацию: для ввода с клавиатуры напишите 'клавиатура',"
                          " для ввода из файла напишите 'файл'\n", ["файл", "клавиатура"])
u2 = typeofread(False, "", "Выберете куда выводить информацию: для вывода на клавиатуру напишите 'клавиатура',"
                           " для вывода в файл напишите 'файл'\n", ["файл", "клавиатура"])
fi = False
if (u == "файл"):
    fi = True
    file = input("Название файла из которого считываются данные\n")
    while not os.path.exists(file):
        print("Файла не существует")
        file = input("Название файла из которого считываются данные\n")
f2 = ""
if (u2 == "файл"):
    fi2 = True
    file2 = input("Название файла куда записать данные\n")
else:
    fi2 = False
if fi:
    fil = open(file, 'r')
if fi2:
    f2 = open(file2, "w")
n = int(input("Введите количество точек(от 8 до 12):\n"))
x1 = [0] * n
y1 = [0] * n
xy = [0] * n
x2 = [0] * n
y2 = [0] * n
x2y = [0] * n
x3y = [0] * n
lnx, lny, lnxy, ylnx, xlny = [0] * n, [0] * n, [0] * n, [0] * n, [0] * n
if fi:
    x1 = list(map(float, fil.readline().split(" ")))
    y1 = list(map(float, fil.readline().split(" ")))
else:
    print("Введите x:\n")
    x1 = list(map(float, input().split()))
    print("Введите y:\n")
    y1 = list(map(float, input().split()))
for i in range(n):
    xy[i] = x1[i] * y1[i]
    x2[i] = x1[i] ** 2
    y2[i] = y1[i] ** 2
    x2y[i] = x2[i] * y1[i]
    x3y[i] = x2y[i] * x1[i]
    lnx[i] = ln(x1[i])
    lny[i] = ln(y1[i])
    lnxy[i] = lnx[i] * lny[i]
    ylnx[i] = y1[i] * lnx[i]
    xlny[i] = x1[i] * lny[i]
sx1 = sum(x1) / n
sy1 = sum(y1) / n
sz, scx, scy = 0, 0, 0
for i in range(n):
    sz += (x1[i] - sx1) * (y1[i] - sy1)
    scx += (x1[i] - sx1) ** 2
    scy += (y1[i] - sy1) ** 2
r = sz / sqrt(scx * scy)
R2 = [0] * 6
fr = [0] * 6
a1 = round((sum(x1) * sum(y1) - n * sum(xy)) / (sum(x1) ** 2 - n * sum(x2)), 4)
b1 = round((sum(x1) * sum(xy) - sum(x2) * sum(y1)) / (sum(x1) ** 2 - n * sum(x2)), 4)
fr[0] = a1 * x + b1
printf("Коэффициент корреляции Пирсона: " + str(r), fi2, f2)
R2 = fun(fr[0], x1, y1, [a1, b1], fi2, f2, R2,0)
Matrix = [[sum(x2), sum(x1), n], [sum(map(lambda x: x ** 3, x1)), sum(x2), sum(x1)],
          [sum(map(lambda x: x ** 4, x1)), sum(map(lambda x: x ** 3, x1)), sum(x2)]]
B = [sum(y1), sum(xy), sum(x2y)]
ab = M(Matrix, B)
fr[1] = ab[0] * x ** 2 + ab[1] * x + ab[2]
R2 = fun(fr[1], x1, y1, ab, fi2, f2, R2, 1)
Matrix = [[sum(map(lambda x: x ** 3, x1)), sum(x2), sum(x1), n],
          [sum(map(lambda x: x ** 4, x1)), sum(map(lambda x: x ** 3, x1)), sum(x2), sum(x1)]
    , [sum(map(lambda x: x ** 5, x1)), sum(map(lambda x: x ** 4, x1)), sum(map(lambda x: x ** 3, x1)), sum(x2)],
          [sum(map(lambda x: x ** 6, x1)), sum(map(lambda x: x ** 5, x1)), sum(map(lambda x: x ** 4, x1)),
           sum(map(lambda x: x ** 3, x1))]]
B = [sum(y1), sum(xy), sum(x2y), sum(x3y)]
abc = M(Matrix, B)
fr[2] = abc[0] * x ** 3 + abc[1] * x ** 2 + abc[2] * x + abc[3]
R2 = fun(fr[2], x1, y1, abc, fi2, f2, R2, 2)
if min(x1) > 0:
    a = round((n * sum(xlny) - sum(x1) * sum(lny)) / (n * sum(map(lambda x: x ** 2, x1)) - sum(x1) ** 2), 4)
    b = round((1 / n * sum(lny) - a / n * sum(x1)), 4)
    fr[3] = exp(a * x + b)
    R2 = fun(fr[3], x1, y1, [a, b], fi2, f2, R2, 3)
    b = round((n * sum(ylnx) - sum(lnx) * sum(y1)) / (n * sum(map(lambda x: x ** 2, lnx)) - sum(lnx) ** 2), 4)
    a = round(1 / n * sum(y1) - b / n * sum(lnx), 4)
    fr[4] = a + b * ln(x)
    R2 = fun(fr[4], x1, y1, [a, b], fi2, f2, R2, 4)
    b = round((n * sum(lnxy) - sum(lnx) * sum(lny)) / (n * sum(map(lambda x: x ** 2, lnx)) - sum(lnx) ** 2), 4)
    a = round(exp(1 / n * sum(lny) - b / n * sum(lnx)), 4)
    fr[5] = a * x ** b
    fun(fr[5], x1, y1, [a, b], fi2, f2, R2, 5)
for i in range(6):
    if R2[i] == max(R2):
        printf("Наилучшая аппроксимирующая функция:" + str(fr[i]), fi2, f2)
        break


def graf(x1, y1, fr, j):
    if j == 5:
        mn = max(min(x1) - 2, 0.001)
        mx = max(x1) + 2
        x213 = np.arange(mn, mx, 0.1)
    else:
        x213 = np.arange(min(x1)*0.9, max(x1)*1.1, 0.1)
    frc = lambdify(x, str(fr))
    if min(x1) > 0:
        plt.subplot(2, 3, j)
    else:
        plt.subplot(1, 3, j)
    plt.plot(x213, frc(x213))
    plt.scatter(x1, y1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'y={fr}')
    plt.grid(True)


if min(x1) > 0:
    fig, axs = plt.subplots(2, 3)
else:
    fig, axs = plt.subplots(1, 3, figsize=(12, 12))
graf(x1, y1, fr[0], 1)
graf(x1, y1, fr[1], 2)
graf(x1, y1, fr[2], 3)
if min(x1) > 0:
    graf(x1, y1, fr[3], 4)
    graf(x1, y1, fr[4], 5)
    graf(x1, y1, fr[5], 6)
plt.show()