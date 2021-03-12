package linearSystems

import (
	"math"
)

const kMax = 100
const eps1 = 0.01
const eps2 = 0.00001
const eps = eps2

func JacobiMethod(A [][]float64, b []float64) ([]float64, int) {
	// вектор решения
	x := make([]float64, len(b))

	// задать начальное приближение
	for i := range x {
		x[i] = 1
	}
	// промежуточный вектор решения
	tx := make([]float64, len(x))

	var k int
	for k = 1; findVectorNorma(subtractVectors(x, tx)) > eps; k++ {
		copy(tx, x)
		for i := range A {
			sum := 0.0
			for j := range A[i] {
				if j != i {
					sum += A[i][j] * tx[j]
				}
			}
			x[i] = (b[i] - sum) / A[i][i]
		}
	}
	return x, k
}

func GaussZeidelMethod(A [][]float64, b []float64) ([]float64, int) {
	// вектор решения
	x := make([]float64, len(b))

	// задать начальное приближение
	for i := range x {
		x[i] = 1
	}
	// промежуточный вектор решения
	tx := make([]float64, len(x))

	var k int
	for k = 1; findVectorNorma(subtractVectors(x, tx)) > eps; k++ {
		copy(tx, x)
		for i := range A {
			sum := 0.0
			for j := range A[i] {
				if j != i {
					sum += A[i][j] * x[j]
				}
			}
			x[i] = (b[i] - sum) / A[i][i]
		}
	}
	return x, k
}

func SteepestDescentMethod(A [][]float64, b []float64) ([]float64, int) {
	// вектор решения
	x := make([]float64, len(b))

	// задать начальное приближение
	for i := range x {
		x[i] = 1
	}
	var k int
	for k = 1; k <= kMax; k++ {
		r := vectorsSubtract(b, matrixVectorMultiply(A, x))
		alpha := dot(r, r) / dot(matrixVectorMultiply(A, r), r)
		t := vectorsSubtract(matrixVectorMultiply(A, x), b)
		x = vectorsSubtract(x, numberVectorMultiply(t, alpha))
		if alpha * findVectorNorma(r) < eps {
			break
		}
	}
	return x, k
}

func MinResidualMethod(A [][]float64, b []float64) ([]float64, int) {
	// вектор решения
	x := make([]float64, len(b))

	// задать начальное приближение
	for i := range x {
		x[i] = 1
	}
	// вектор невязки
	r := make([]float64, len(b))
	copy(r, x)

	var k int
	for k = 1; k <= kMax && findVectorNorma(r) > eps; k++ {
		r = vectorsSubtract(b, matrixVectorMultiply(A, x))
		t := matrixVectorMultiply(A, r)
		tau := dot(t, r) / (findVectorNorma(t) * findVectorNorma(t))
		x = vectorsSum(x, numberVectorMultiply(r, tau))
	}
	return x, k
}

func CheckConvergenceCondition(A [][]float64) bool {
	for i := range A {
		sum := 0.0
		for j := range A[i] {
			if i != j {
				sum += math.Abs(A[i][j])
			}
		}
		if math.Abs(A[i][i]) <= sum {
			return false
		}
	}
	return true
}

func findVectorNorma(v []float64) float64 {
	sum := 0.0
	for i := range v {
		sum += v[i] * v[i]
	}
	return math.Sqrt(sum)
}

func subtractVectors(v []float64, w []float64) []float64 {
	res := make([]float64, len(v))
	for i := range v {
		res[i] = v[i] - w[i]
	}
	return res
}

func matrixVectorMultiply(A [][]float64, b []float64) []float64 {
	v := make([]float64, len(b))
	for i := range A {
		for j := range A[i] {
			v[i] += A[i][j] * b[j]
		}
	}
	return v
}

func dot(a []float64, b []float64) float64 {
	x := 0.0
	for i := range a {
		x += a[i] * b[i]
	}
	return x
}

func vectorsSum(a []float64, b []float64) []float64 {
	v := make([]float64, len(a))
	for i := range a {
		v[i] = a[i] + b[i]
	}
	return v
}

func vectorsSubtract(a []float64, b []float64) []float64 {
	v := make([]float64, len(a))
	for i := range a {
		v[i] = a[i] - b[i]
	}
	return v
}

func numberVectorMultiply(x []float64, num float64) []float64 {
	v := make([]float64, len(x))
	for i := range x {
		v[i] = x[i] * num
	}
	return v
}