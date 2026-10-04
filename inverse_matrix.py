"""
역행렬 계산 프로그램 (이산수학 과제 1)

 - 방법 1: 행렬식(여인수 전개 + 수반행렬)을 이용한 역행렬 계산
 - 방법 2: 가우스-조던 소거법을 이용한 역행렬 계산
 - 두 방법의 결과를 출력하고 동일한지 비교

모든 계산은 fractions.Fraction(유리수)으로 수행하여 부동소수점 오차 없이
정확한 값을 구한다.
"""
from fractions import Fraction
import random
import time


class SingularMatrixError(Exception):
    """역행렬이 존재하지 않는(특이) 행렬일 때 발생시키는 예외"""


# ============================================================
# 공통 유틸리티
# ============================================================
def fmt(x):
    """Fraction을 보기 좋은 문자열로 변환 (정수면 정수, 아니면 a/b)"""
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def print_matrix(m, title=None, indent="  ", decimal=False):
    """행렬을 열 너비를 맞추어 출력"""
    if title:
        print(f"{indent}{title}")
    cells = [[f"{float(x):.4f}" if decimal else fmt(x) for x in row] for row in m]
    width = max(len(c) for row in cells for c in row)
    for row in cells:
        print(f"{indent}[ " + "  ".join(c.rjust(width) for c in row) + " ]")


def print_augmented(aug, n, title):
    """첨가행렬 [A | I]를 가운데 구분선과 함께 출력"""
    print(f"    {title}")
    cells = [[fmt(x) for x in row] for row in aug]
    width = max(len(c) for row in cells for c in row)
    for row in cells:
        left = "  ".join(c.rjust(width) for c in row[:n])
        right = "  ".join(c.rjust(width) for c in row[n:])
        print(f"    [ {left} | {right} ]")


def identity(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def transpose(m):
    return [list(col) for col in zip(*m)]


def multiply(a, b):
    """행렬 곱 a × b"""
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


# ============================================================
# 1. 행렬 입력 기능
# ============================================================
def input_size():
    """정방행렬의 크기 n(1 이상의 정수)을 입력받는다."""
    while True:
        s = input("행렬의 크기 n을 입력하세요: ").strip()
        try:
            n = int(s)
            if n >= 1:
                return n
        except ValueError:
            pass
        print("  [입력 오류] 1 이상의 정수를 입력하세요.")


def input_matrix(n):
    """n × n 행렬을 한 행씩 입력받아 2차원 리스트로 반환한다."""
    print(f"{n} x {n} 행렬을 한 행씩 입력하세요. (공백 구분 / 정수·소수·분수 가능 예: 3 -1 0.5 2/3)")
    matrix = []
    while len(matrix) < n:
        tokens = input(f"  {len(matrix) + 1}행: ").split()
        if len(tokens) != n:
            print(f"  [입력 오류] {n}개의 수를 입력해야 합니다. (입력한 개수: {len(tokens)}개)")
            continue
        try:
            row = [Fraction(t) for t in tokens]
        except (ValueError, ZeroDivisionError):
            print("  [입력 오류] 숫자가 아닌 값이 있습니다. 해당 행을 다시 입력하세요.")
            continue
        matrix.append(row)
    return matrix


def random_matrix(n, low=-9, high=9):
    """[추가 기능] -9 ~ 9 범위의 정수로 이루어진 n × n 랜덤 행렬 생성"""
    return [[Fraction(random.randint(low, high)) for _ in range(n)] for _ in range(n)]


# ============================================================
# 2. 행렬식을 이용한 역행렬 계산 기능
# ============================================================
def get_minor(m, row, col):
    """row행과 col열을 제거한 소행렬(minor matrix)"""
    return [r[:col] + r[col + 1:] for k, r in enumerate(m) if k != row]


def determinant(m):
    """여인수 전개(라플라스 전개)로 행렬식을 재귀적으로 계산"""
    n = len(m)
    if n == 1:
        return m[0][0]
    if n == 2:
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]
    det = Fraction(0)
    for j in range(n):                      # 1행 기준 전개
        if m[0][j] != 0:
            det += (-1) ** j * m[0][j] * determinant(get_minor(m, 0, j))
    return det


def inverse_by_determinant(m):
    """A^-1 = adj(A) / det(A)  (adj(A): 여인수행렬의 전치행렬)

    반환값: (역행렬, 행렬식, 수반행렬)
    """
    n = len(m)
    det = determinant(m)
    if det == 0:
        raise SingularMatrixError("행렬식 det(A) = 0 이므로 역행렬이 존재하지 않습니다.")
    if n == 1:
        adj = [[Fraction(1)]]
    else:
        cofactor = [[(-1) ** (i + j) * determinant(get_minor(m, i, j))
                     for j in range(n)] for i in range(n)]
        adj = transpose(cofactor)
    inv = [[adj[i][j] / det for j in range(n)] for i in range(n)]
    return inv, det, adj


# ============================================================
# 3. 가우스-조던 소거법을 이용한 역행렬 계산 기능
# ============================================================
def inverse_by_gauss_jordan(m, show_steps=False):
    """첨가행렬 [A | I]에 기본 행 연산을 적용하여 [I | A^-1]로 만든다."""
    n = len(m)
    aug = [list(m[i]) + identity(n)[i] for i in range(n)]
    if show_steps:
        print_augmented(aug, n, "초기 첨가행렬 [A | I]")

    for col in range(n):
        # 부분 피벗팅: 현재 열에서 절댓값이 가장 큰 원소를 가진 행을 피벗으로 선택
        pivot_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if aug[pivot_row][col] == 0:
            raise SingularMatrixError(
                f"{col + 1}열에서 0이 아닌 피벗을 찾을 수 없습니다. 역행렬이 존재하지 않습니다.")
        ops = []
        if pivot_row != col:                # (1) 행 교환
            aug[col], aug[pivot_row] = aug[pivot_row], aug[col]
            ops.append(f"R{col + 1} <-> R{pivot_row + 1}")
        pivot = aug[col][col]
        if pivot != 1:                      # (2) 피벗을 1로 만들기
            aug[col] = [x / pivot for x in aug[col]]
            ops.append(f"R{col + 1} <- R{col + 1} / ({fmt(pivot)})")
        for r in range(n):                  # (3) 피벗 열의 나머지 원소를 0으로
            if r != col and aug[r][col] != 0:
                factor = aug[r][col]
                aug[r] = [a - factor * b for a, b in zip(aug[r], aug[col])]
                ops.append(f"R{r + 1} <- R{r + 1} - ({fmt(factor)})R{col + 1}")
        if show_steps:
            print(f"    ({col + 1}단계) " + ", ".join(ops))
            print_augmented(aug, n, f"{col + 1}열 소거 후")

    return [row[n:] for row in aug]


# ============================================================
# 4. 결과 출력 및 비교 기능
# ============================================================
def matrices_equal(a, b):
    return all(x == y for ra, rb in zip(a, b) for x, y in zip(ra, rb))


def compare_results(inv_det, inv_gj):
    print("\n[4] 결과 비교")
    if inv_det is None and inv_gj is None:
        print("  => 두 방법 모두 '역행렬이 존재하지 않음'으로 판정하였습니다. (결과 일치)")
    elif inv_det is None or inv_gj is None:
        print("  => 한 방법만 역행렬을 구했습니다. 결과가 일치하지 않습니다!")
    elif matrices_equal(inv_det, inv_gj):
        print("  => 두 방법으로 구한 역행렬이 모든 원소에서 동일합니다. (결과 일치)")
    else:
        diff = max(abs(x - y) for ra, rb in zip(inv_det, inv_gj) for x, y in zip(ra, rb))
        print(f"  => 두 결과가 다릅니다. (최대 원소 차이: {float(diff):.6g})")


def verify_inverse(m, inv):
    """[추가 기능] A × A^-1 이 단위행렬인지 검증"""
    product = multiply(m, inv)
    print("\n[5] 검증 (추가 기능): A x A^-1 = I ?")
    print_matrix(product, "A x A^-1 =")
    if matrices_equal(product, identity(len(m))):
        print("  => 단위행렬 I 가 나오므로 올바른 역행렬입니다.")
    else:
        print("  => 단위행렬이 아닙니다. 계산 오류!")


# ============================================================
# 전체 실행 흐름
# ============================================================
def solve(m, show_steps):
    n = len(m)
    print()
    print_matrix(m, f"입력 행렬 A ({n} x {n})", indent="")

    # 방법 1: 행렬식
    print("\n[2] 행렬식을 이용한 역행렬")
    inv_det = None
    t0 = time.perf_counter()
    try:
        inv_det, det, adj = inverse_by_determinant(m)
        t_det = time.perf_counter() - t0
        print(f"  det(A) = {fmt(det)}")
        print_matrix(adj, "수반행렬 adj(A) =")
        print_matrix(inv_det, "A^-1 = adj(A) / det(A) =")
        print(f"  (계산 시간: {t_det * 1000:.3f} ms)")
    except SingularMatrixError as e:
        print(f"  [오류] {e}")

    # 방법 2: 가우스-조던
    print("\n[3] 가우스-조던 소거법을 이용한 역행렬")
    inv_gj = None
    t0 = time.perf_counter()
    try:
        inv_gj = inverse_by_gauss_jordan(m, show_steps)
        t_gj = time.perf_counter() - t0
        print_matrix(inv_gj, "A^-1 =")
        if not show_steps:
            print(f"  (계산 시간: {t_gj * 1000:.3f} ms)")
    except SingularMatrixError as e:
        print(f"  [오류] {e}")

    compare_results(inv_det, inv_gj)
    if inv_gj is not None:
        verify_inverse(m, inv_gj)
        print()
        print_matrix(inv_gj, "[참고] A^-1 소수 표현:", decimal=True)


def benchmark():
    """[추가 기능] 랜덤 행렬로 두 방법의 수행 시간 비교"""
    print("\n[추가 기능] n = 2 ~ 7 랜덤 행렬에 대한 계산 시간 비교")
    print("    n |   행렬식(ms) | 가우스-조던(ms) | 결과 일치")
    print("  ----+--------------+-----------------+----------")
    for n in range(2, 8):
        m = random_matrix(n)
        while determinant(m) == 0:
            m = random_matrix(n)
        t0 = time.perf_counter()
        a, _, _ = inverse_by_determinant(m)
        t1 = time.perf_counter()
        b = inverse_by_gauss_jordan(m)
        t2 = time.perf_counter()
        same = "O" if matrices_equal(a, b) else "X"
        print(f"  {n:>3} | {(t1 - t0) * 1000:>12.3f} | {(t2 - t1) * 1000:>15.3f} | {same:>5}")


def main():
    while True:
        print("=" * 50)
        print("        역행렬 계산 프로그램 (Inverse Matrix)")
        print("=" * 50)
        print(" 1. 행렬 직접 입력하여 역행렬 구하기")
        print(" 2. 랜덤 행렬로 역행렬 구하기 (추가 기능)")
        print(" 3. 두 방법의 계산 시간 비교 (추가 기능)")
        print(" 0. 종료")
        choice = input("메뉴 선택 >> ").strip()

        if choice == "0":
            print("프로그램을 종료합니다.")
            break
        elif choice == "1":
            print("\n[1] 행렬 입력")
            n = input_size()
            m = input_matrix(n)
        elif choice == "2":
            print("\n[1] 랜덤 행렬 생성")
            n = input_size()
            m = random_matrix(n)
        elif choice == "3":
            benchmark()
            print()
            continue
        else:
            print("  [입력 오류] 0~3 중에서 선택하세요.\n")
            continue

        steps = input("가우스-조던 소거 과정을 단계별로 출력할까요? (y/n): ").strip().lower()
        solve(m, steps == "y")
        print()


if __name__ == "__main__":
    main()
