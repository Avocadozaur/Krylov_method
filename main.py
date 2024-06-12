import math

import numpy as np


def Kuchta_Wiktor_(M):
    def Gaussior(A, b):
        n = len(b)
        for i in range(n):
            for j in range(i + 1, n):
                factor = A[j, i] / A[i, i]
                for k in range(i, n):
                    A[j, k] -= factor * A[i, k]
                b[j] -= factor * b[i]
        # Tu mam macierz gornotrojkatna
        x = np.zeros(n)
        # od dolu obliczam wartosci
        for i in range(n - 1, -1, -1):
            x[i] = b[i]
            for j in range(i + 1, n):
                x[i] -= A[i, j] * x[j]
            x[i] /= A[i, i]
        return x

    b = np.array([1, 0])
    for i in range(len(M[0]) - 2):
        b = np.insert(b, 1, 0)  # dostosowanie wektora b do wielkosci macierzy
    Ab = np.matmul(b, M)
    A2b = np.matmul(Ab, M)
    A3b = np.matmul(A2b, M)
    A = np.matrix([np.ravel(A2b), np.ravel(Ab), b])
    coeffs = Gaussior(np.transpose(A), -np.ravel(A3b))
    coeffs = np.insert(coeffs, 0, 1)  # wrzucam jeden na początek zgodnie ze wzorem

    def F_X(lista_coeffs, iks):
        wynik = 0
        for i in range(len(lista_coeffs)):
            wynik += lista_coeffs[i] * math.pow(iks, len(lista_coeffs) - 1 - i)
        return wynik

    def sieczne(x0, x1, interations):
        x2 = 0
        for i in range(interations):
            if abs(x1 - x0) < 0.0001:  # dokładnosc, np 1/10 000
                return x2
            x2 = ((F_X(coeffs, x1) * x0 - F_X(coeffs, x0) * x1) / (
                    F_X(coeffs, x1) - F_X(coeffs, x0)))
            x0, x1 = x1, x2
        return x2

    def horner(lista_coeffs, root):
        n = len(lista_coeffs)
        iloraz = [0] * (n - 1)
        reszta = lista_coeffs[0]

        for i in range(1, n):
            iloraz[i - 1] = reszta
            reszta = lista_coeffs[i] + reszta * root

        return iloraz

    lista = []
    while len(coeffs) >= 2:
        lista.append(sieczne(0, 1, 10))
        coeffs = horner(coeffs, sieczne(0, 1, 10))

    wynik = ''.join([f"(x - {pierwiastek}) " for pierwiastek in lista])
    return wynik


M = [[2.0, 1.0, 2.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]]
# M = [[1., -2., 3.], [-2., 2., -1.], [3., -1., 0.]]
print(Kuchta_Wiktor_(M))
M = [[2.0, 1.0], [1.0, 2.0]]
print(Kuchta_Wiktor_(M))
# M = [[2.0, 1.0, 0.0, 3.0], [1.0, 4.0, 5.0, 1.0], [0.0, 5.0, 6.0, 2.0], [3.0, 1.0, 2.0, 7.0]]
# print(Kuchta_Wiktor_(M))
