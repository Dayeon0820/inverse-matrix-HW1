# 역행렬 계산 프로그램 (Inverse Matrix Calculator)

이산수학 과제 1 — n × n 정방행렬의 역행렬을 **행렬식(수반행렬)** 과 **가우스-조던 소거법** 두 가지 방법으로 계산하고 결과를 비교하는 Python 프로그램입니다.

## 실행 방법

```bash
git clone https://github.com/Dayeon0820/inverse-matrix-HW1.git
cd inverse-matrix-HW1
python inverse_matrix.py
```

- Python 3.8 이상, 표준 라이브러리만 사용 (`fractions`, `random`, `time`)

## 기능

| 구분 | 기능 | 함수 |
|---|---|---|
| 구성 조건 1 | n 입력 후 n × n 행렬을 행 단위로 입력, 2차원 리스트로 저장 | `input_size`, `input_matrix` |
| 구성 조건 2 | 행렬식(여인수 전개)과 수반행렬로 역행렬 계산, det = 0이면 오류 메시지 | `determinant`, `inverse_by_determinant` |
| 구성 조건 3 | 가우스-조던 소거법으로 역행렬 계산, 피벗이 0이면 예외 처리 | `inverse_by_gauss_jordan` |
| 구성 조건 4 | 두 결과 출력 및 동일 여부 비교 | `compare_results` |
| 추가 기능 | 분수·소수 입력과 정확한 유리수(Fraction) 연산 | `input_matrix`, `fmt` |
| 추가 기능 | 잘못된 입력(크기·개수·문자) 검증 후 재입력 | `input_size`, `input_matrix` |
| 추가 기능 | 가우스-조던 행 연산 과정 단계별 출력 | `inverse_by_gauss_jordan(show_steps=True)` |
| 추가 기능 | A × A⁻¹ = I 검증, 소수 표현 출력 | `verify_inverse` |
| 추가 기능 | 랜덤 행렬 생성, 두 방법의 계산 시간 비교 | `random_matrix`, `benchmark` |

## 실행 예시

```
[1] 행렬 입력
행렬의 크기 n을 입력하세요: 2
  1행: 4 7
  2행: 2 6
...
[2] 행렬식을 이용한 역행렬
  det(A) = 10
  A^-1 = adj(A) / det(A) =
  [   3/5  -7/10 ]
  [  -1/5    2/5 ]
[4] 결과 비교
  => 두 방법으로 구한 역행렬이 모든 원소에서 동일합니다. (결과 일치)
```
