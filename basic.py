import numpy
import scipy.integrate as scipy

N = int(input("dawaj dokladnosc przedzialu kjurwa"))
u_2 = 0
domain = [0, 2]


def E(x):
    if x <= 1 and x >= 0:
        return 3
    elif x > 1 and x <= 2:
        return 5
    else:
        return 0


def solveEquationsMatrixes(B, L):
    print(numpy.array(B))
    print(numpy.array(L))
    return numpy.linalg.solve(B, L)


def getX(n):
    return 2*n/N


def e(x, i):
    if i > 0:
        if x < getX(i - 1) or x > getX(i + 1):
            return 0
        if x < getX(i):
            return N/2 * (x - getX(i-1))
        else:
            return N/2 * (getX(i+1) - x)
    elif i == 0:
        if x > getX(i + 1):
            return 0
        else:
            return N/2 * (getX(i+1) - x)


def eFun(i):
    return lambda x: e(x, i)


def diffe(x, i):
    if i > 0:
        if x <= getX(i-1) or x >= getX(i+1):
            return 0
        elif x < getX(i):
            return N/2
        else:
            return -N/2
    elif i == 0:
        if x > getX(i+1):
            return 0
        else:
            return -N/2

def diffeFun(i):
    return lambda x: diffe(x, i)

def integrate(f):
    return 3*scipy.quad(f,0,1)[0] + 5*scipy.quad(f,1,2)[0]


# def integral(i, j):
#     # if abs(i - j) > 1:
#     #     return 0
#     # elif i != j:
#     #     xl = getX(min(i, j))
#     #     xr = xl + 2/N
#     #     if xl > 1:
#     #         return -N*N*5 * (xr - xl)
#     #     elif xr < 1:
#     #         return -N*N*3 * (xr - xl)
#     #     else:
#     #         return -N*N*(3 * (1 - xl) + 5 * (xr - 1))
#     # else:
#     #     xl = getX(i)
#     #     xr = xl + 2/N
#     #     if xl > 1:
#     #         return -N * N * 5 * (xr - xl)
#     #     elif xr < 1:
#     #         return -N * N * 3 * (xr - xl)
#     #     else:
#     #         return -N * N * (3 * (1 - xl) + 5 * (xr - 1))
#      return scipy.quad(lambda x: E(x)*diffe(x, i)*diffe(x, j), 0, 2)[0]



def BLinear(i, j):
    constant = 0
    if i == j and i == 0:
        constant = 1
    return constant + integrate(
        lambda x: diffe(x, i) * diffe(x, j))


def solve():
    LMatrix = []
    BMatrix = numpy.zeros((N+1, N+1))
    for i in range(0, N):
        array = []
        if i == 0:
            LMatrix.append(-20)
        else:
            LMatrix.append(0)
        for j in range(0, N):
            array.append(BLinear(j, i))
        array.append(0)
        BMatrix[i] = array
    LMatrix.append(0)
    array = numpy.zeros(N+1)
    array[N] = 1
    BMatrix[N] = array
    BMatrix = numpy.array(BMatrix)
    return solveEquationsMatrixes(BMatrix, LMatrix)


def getBparameters(i, j):
    return BLinear(j, i)


def getBMatrix():
    BMatrix = numpy.zeros((N+1, N+1))
    for i in range(N+1):
        for j in range(N+1):
            BMatrix[i][j] = getBparameters(i, j)
            if i == N or j == N:
                BMatrix[i][j] = 0
            BMatrix[N][N] = 1
    return BMatrix


def getLMatrix():
    matrix = numpy.zeros((N+1))
    matrix[0] = -30*e(0, 0)
    return matrix

def solveForU():
    factors = solve()
    u = lambda x: sum([factors[i] * e(x, i) for i in range(0, N + 1)])
    return u
