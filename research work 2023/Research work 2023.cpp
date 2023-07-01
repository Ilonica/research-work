#include <iostream>
#include <iomanip>
#include <math.h>
#include <omp.h>
#include <corecrt_math_defines.h>
using namespace std;

#define THREADS_COUNTS omp_get_num_threads()
#define CPU_COUNTS omp_get_num_procs()
#define CPU_ID omp_get_thread_num()
using namespace std;

enum verbose { off, medium, full };
verbose DetailedInf = medium;

const double eps = 1e-5;
//double n0 = 100;
const double π = M_PI;
const double θ = 1.0 / 3.0;
//
//const double a = 0;
//const double b = 5;
//const double X0 = 0.1;
//double f(double x) {
//	return 2 * x + 1;
//}

double n0 = 58715500;
const double a = 0;
const double b = 2;
const double X0 = 3;
double f(double x) {
	return sin(x);
}

//const double a = -1;
//const double b = 10;
//const double X0 = 0.7;
//double f(double x) {
//	return pow(x, 2) + 3 * x - 2;
//}
//
//const double a = 0;
//const double b = 2;
//const double X0 = 0.1;
//double f(double x) {
//	return (1 / pow(M_E, 2 * x));
//}

//const double a = 0;
//const double b = π;
//const double X0 = 0.5;
//double f(double x) {
//	return 4 / (x * x + 1);
//}

double In; double In1; double I2n;

double RectangleMethod1(double b, int n) {
	double sum = 0; double h = (b - a) / n;
	int i;

#pragma omp parallel for reduction(+:sum)
	for (i = 0; i < n; i++) {
		if (DetailedInf == full) {
			cout << "Потоки: " << THREADS_COUNTS << endl;
			cout << CPU_ID << " - " << sum << endl;
		}
		sum += f(a + i * h);
	}
	sum *= h;
	if (DetailedInf) cout << "b = " << fixed << b << " Sum = " << sum << " N = " << n << endl;
	return sum;
}

double RungeRule(double b) {
	double In = RectangleMethod1(b, n0); double In2 = RectangleMethod1(b, 2 * n0);
	double diff = In2 - In;

	while (θ * abs(diff) >= eps) {
		n0 *= 2;
		In = RectangleMethod1(b, n0); In2 = RectangleMethod1(b, 2 * n0);
		diff = In2 - In;
	}
	if (DetailedInf) cout << "N1 = " << n0 << " Diff = " << diff << endl;
	return n0;
}

double ChordMethod() {
	double start = omp_get_wtime();
	double x0 = X0;
	double x1 = b;
	double N0 = RungeRule(x0);
	In = RectangleMethod1(x0, N0);
	In1 = RectangleMethod1(x1, N0);
	double x2 = x1 - ((In1 - b) * (x1 - x0)) / ((In1 - b) - (In - b));
	int iter = 1;

	while (abs(x2 - x1) >= eps) {
		iter++;
		x0 = x1;
		x1 = x2;
		In = In1;
		In1 = RectangleMethod1(x1, N0);
		x2 = x1 - ((In1 - b) * (x1 - x0)) / ((In1 - b) - (In - b));
	}
	double end = omp_get_wtime();
	cout << "Результат: x1 = " << fixed << x0 << endl;
	cout << "Количество итераций: " << iter << endl;

	return end - start;
}

void iData() {
	cout.precision(10);
	cout << "Входные данные:" << endl;
	cout << "Нижний предел интегрирования a = " << a << endl;
	cout << "Верхний предел интегрирования b = " << b << endl;
	cout << "Допустимая погрешность решения eps = " << eps << endl;
	cout << "Количество отрезков n0 = " << n0 << endl;
	cout << "Начальное приближение x0 = " << X0 << endl;
	cout << "Общее количество процессоров: " << CPU_COUNTS << endl;
}

void setCPUcounts(int N) {
	omp_set_num_threads(N);
	cout << "Количество потоков: " << N << endl;
}

void solve() {
	cout << endl;
	cout.precision(10);
	cout << "Время вычисления: " << ChordMethod() << " секунд" << endl;
}

int main() {
	setlocale(0, "");
	iData();
	setCPUcounts(9); // Количество потоков
	solve();
	return 0;
}