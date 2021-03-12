package main

import (
	"fmt"
	"main/LinearSystems"
)

func testJacobiMethod(A [][]float64, b []float64) {
	x, k := linearSystems.JacobiMethod(A, b)
	fmt.Println("J\t", k, "iter.")
	fmt.Println(x)
}

func testGaussZeidelMethod(A [][]float64, b []float64) {
	x, k := linearSystems.GaussZeidelMethod(A, b)
	fmt.Println("GZ\t", k, "iter.")
	fmt.Println(x)
}

func testSteepestDescentMethod(A [][]float64, b []float64) {
	x, k := linearSystems.SteepestDescentMethod(A, b)
	fmt.Println("SD\t", k, "iter.")
	fmt.Println(x)
}

func testMinResidualMethod(A [][]float64, b []float64) {
	x, k := linearSystems.MinResidualMethod(A, b)
	fmt.Println("MR\t", k, "iter.")
	fmt.Println(x)
}

func main() {
	// initial equation
	A := [][]float64{
		{0.341, -0.542, 0.418, -2.110},
		{-0.111, 0.915, 0.012, 0.341},
		{0.546, 0.211, -0.318, 1.810},
		{0.302, 0.201, 2.130, -0.115},
	}
	b := []float64{-1.893, 2.249, 2.249,  2.518}

	isConvergent := linearSystems.CheckConvergenceCondition(A)
	fmt.Println("Convergent:", isConvergent)

	// transformed equation
	A = [][]float64{
		{0.341+0.546, -0.542+0.211, 0.418-0.318, -2.110+1.810},
		{-0.111, 0.915, 0.012, 0.341},
		{0.302, 0.201, 2.130, -0.115},
		{0.546, 0.211, -0.318, 1.810},
	}
	b = []float64{-1.893+2.249, 2.249, 2.518, 2.249}

	isConvergent = linearSystems.CheckConvergenceCondition(A)
	fmt.Println("Convergent:", isConvergent)

	testJacobiMethod(A, b)
	testGaussZeidelMethod(A, b)
	testSteepestDescentMethod(A, b)
	testMinResidualMethod(A, b)

	actualRoots := []float64{
		1.42439450888372,
		2.368335563925662,
		0.7932173571765274,
		0.6761353114603924}
	fmt.Print("roots\n", actualRoots)

	// 1.42439450888372, 2.368335563925662, 0.7932173571765274, 0.6761353114603924
}

