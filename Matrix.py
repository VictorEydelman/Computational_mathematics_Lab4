def M(Matrix,B):
    n = len(Matrix)
    for i in range(n):
        max_el = abs(Matrix[i][i])
        max_row = i
        for k in range(i + 1, n):
            if abs(Matrix[k][i]) > max_el:
                max_el = abs(Matrix[k][i])
                max_row = k

        Matrix[i], Matrix[max_row] = Matrix[max_row], Matrix[i]
        B[i], B[max_row] = B[max_row], B[i]

        for k in range(i + 1, n):
            c = -Matrix[k][i] / Matrix[i][i]
            for j in range(i, n):
                if i == j:
                    Matrix[k][j] = 0
                else:
                    Matrix[k][j] += c * Matrix[i][j]
            B[k] += c * B[i]

    x = [0 for _ in range(n)]
    for i in range(n - 1, -1, -1):
        x[i] = round(B[i] / Matrix[i][i],4)
        for k in range(i - 1, -1, -1):
            B[k] -= Matrix[k][i] * x[i]
    return x