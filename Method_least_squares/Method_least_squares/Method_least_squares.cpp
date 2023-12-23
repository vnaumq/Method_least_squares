#include <iostream>

using namespace std;

double* Gauss(double** Array_A, double* Array_B, int n);

int main()
{
	setlocale(LC_ALL, "rus");
	const int N = 10;
	const int m = 2;
	double x[10] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
	double y[10] = { 0.957, 0.969, 0.976, 0.978, 0.975, 0.968, 0.954, 0.939, 0.918, 0.894 };
	double* POWERX = new double[2 * m];
	cout << "POWERX: ";
	for (int k = 0; k < 2 * m; k++)
	{
		POWERX[k] = 0;
		for (int i = 0; i < N; i++)
		{
			POWERX[k] += pow(x[i], k + 1);
		}
	}

	double** SUMX = new double* [m + 1];
	for (int i = 0; i < m + 1; i++)
	{
		SUMX[i] = new double[m + 1];
	}
	for (int l = 0; l < m + 1; l++)
	{
		for (int j = 0; j < m + 1; j++)
		{
			if (j + l)
			{
				SUMX[l][j] = POWERX[l + j - 1];
			}
			else
			{
				SUMX[l][j] = N;
			}
		}
	}
	double* PRAW = new double[m + 1];
	for (int l = 0; l < m + 1; l++)
	{
		PRAW[l] = 0;
		for (int i = 0; i < N; i++)
		{
			PRAW[l] += y[i] * pow(x[i], l);
		}
	}
	double* a = Gauss(SUMX, PRAW, m + 1);
	double S2 = 0;
	for (int i = 0; i < N; i++)
	{
		double sum = y[i];
		for (int j = 0; j < m + 1; j++)
		{
			sum -= a[j] * pow(x[i], j);
		}
		S2 += pow(sum, 2);
	}
	S2 /= N - m - 1;
	double sigma = sqrt(S2);
	cout << "Коэфициенты a: " << endl;
	for (int i = 0; i < m + 1; i++)
	{
		cout << a[i] << " ";
	}
	cout << "\nCреднеквадратическое отклонение: " << sigma << endl;
	system("pause");
	return 0;
}

double* Gauss(double** Array_A, double* Array_B, int n)
{
	double* X = new double[n];  //массив ответов
	for (int k = 0; k < n; k++) // прямой ход
	{
		for (int i = k + 1; i < n; i++)
		{
			if (abs(Array_A[i][k]) > abs(Array_A[k][k]))
			{
				swap(Array_A[i], Array_A[k]);  //меняются адреса т.к. двумерный массив
				swap(Array_B[i], Array_B[k]);   //меняются значения
			}
		}
		double A_Main = Array_A[k][k];
		if (A_Main == 0)
		{
			cout << "error\n";
			system("pause");
			exit(0);
		}
		for (int i = k; i < n; i++)
		{
			Array_A[k][i] /= A_Main;
		}
		Array_B[k] /= A_Main;
		for (int i = k + 1; i < n; i++)
		{
			double s = Array_A[i][k];
			for (int j = k; j < n; j++)
			{
				Array_A[i][j] -= s * Array_A[k][j];
			}
			Array_B[i] -= s * Array_B[k];
		}
	}
	for (int k = n - 1; k >= 0; k--)  //обратный ход
	{
		X[k] = Array_B[k];
		for (int i = n - 1; i > k; i--)
		{
			X[k] -= Array_A[k][i] * X[i];
		}
	}
	return X;
}