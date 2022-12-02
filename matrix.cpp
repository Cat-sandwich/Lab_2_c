#pragma once
#include "matrix.h"
#include "exceptions.h"
#include <iostream>
#include <random>
#include <complex>
using namespace std;


Matrix::Matrix()
{
	m = 0;
	n = 0;
	data = NULL;

}

Matrix::Matrix(int m, int n)
{
	this->m = m;
	this->n = n;

	data = (double**) new double* [m];

	for (int i = 0; i < m; i++)
		data[i] = (double*) new double[n];

}


Matrix::Matrix(int m, int n, const double& value)
{
	//*this = Matrix(m, n);
	Set_m(m);
	Set_n(n);

	data = (double**) new double* [m];

	for (int i = 0; i < m; i++)
		data[i] = (double*) new double[n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			data[i][j] = value;
	}

}

Matrix::Matrix(const Matrix& Matrix)
{
	m = Matrix.m;
	n = Matrix.n;


	data = (double**) new double* [m];

	for (int i = 0; i < m; i++)
		data[i] = (double*) new double[n];

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			data[i][j] = Matrix.data[i][j];
}

void Matrix::Set_m(int m)
{
	if (m > 0)
		this->m = m;
}

void Matrix::Set_n(int n)
{
	if (n > 0)
		this->n = n;
}

int Matrix::Get_m()
{
	return m;
}

int Matrix::Get_n()
{
	return n;
}

void Matrix::Set_Data(const double& value)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			data[i][j] = value;
}

double Matrix::Get_Data(int i, int j) const //todo+
{
	if ((m <= i) && (n <= j)) throw Invalid_Index();
	return data[i][j];

}



void Matrix::Print(const int& Number_Matrix)
{
	if (data == NULL) throw Empty();
	cout << "Number_Matrix: " << Number_Matrix << endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			cout << data[i][j] << "\t";
		cout << endl;
	}
}


Matrix& Matrix::operator = (const Matrix& M) //ToDO+
{
	if ((m == M.m) && (n == M.n))
	{
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				data[i][j] = M.data[i][j];
		return *this;
	}
	if ((m < M.m) || (n < M.n))
	{
		for (int i = 0; i < m; i++)
			delete[] data[i];

		delete[] data;

	}
	m = M.m;
	n = M.n;


	data = (double**) new double* [m];

	for (int i = 0; i < m; i++)
		data[i] = (double*) new double[n];

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			data[i][j] = M.data[i][j];
	return *this;
}

Matrix::~Matrix()
{
	if (n > 0)
	{
		for (int i = 0; i < m; i++)
			delete[] data[i];
	}

	if (m > 0)
		delete[] data;
}

double& Matrix::operator ()(int m, int n) const
{
	if ((m > this->m) || (n > this->n)) throw Invalid_Index();

	return   data[m][n];
}

Matrix& Matrix::operator () (int m, int n, const double& value)
{
	if ((m > this->m) || (n > this->n)) throw Invalid_Index();
	data[m][n] = value;
	return *this;

}

Matrix Matrix::operator + (const Matrix& New_Matrix) {
	if (this->m != New_Matrix.m || this->n != New_Matrix.n) throw Different_Dimensions();
	Matrix res(m, n, 0);
	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < n; j++) {
			res.data[i][j] = this->data[i][j] + New_Matrix.data[i][j];
		}
	}
	return res;
}


Matrix Matrix::operator - (const Matrix& New_Matrix) {
	if (this->m != New_Matrix.m || this->n != New_Matrix.n) throw Different_Dimensions();

	Matrix res(m, n, 0);
	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < n; j++) {
			res.data[i][j] = this->data[i][j] - New_Matrix.data[i][j];
		}
	}
	return res;
}


Matrix Matrix::operator * (const Matrix& New_Matrix)
{
	if (n != New_Matrix.m) throw Different_Dimensions();

	Matrix res(m, New_Matrix.n, 0);

	for (int i = 0; i < res.m; ++i)
		for (int j = 0; j < res.n; ++j)
		{

			for (int k = 0; k < m; ++k)

				res.data[i][j] += data[i][k] * New_Matrix.data[k][j];
		}
	return res;
}

Matrix Matrix::operator * (const double& scalar)
{
	Matrix res(m, n, 0);

	for (int i = 0; i < this->m; i++)
	{
		for (int j = 0; j < this->n; j++)
		{
			res.data[i][j] = data[i][j] * scalar;
		}
	}
	return res;
}


Matrix Matrix::operator / (const double& scalar)
{
	if (scalar == (0)) throw Divizion_By_Zero();
	Matrix res(m, n, 0);

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++) {
			res.data[i][j] = data[i][j] / scalar;
		}
	}
	return res;
}

double Matrix::Сalculating_trace_matrix()
{
	if (n != m) throw Dimensions_Incorrect();
	double trace = 0;
	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++) {
			if (i == j)
				trace += data[i][j];
		}
	}
	return trace;
}


Matrix Matrix::Transpose()
{
	Matrix Transposed(m, n);
	for (int i = 0; i < m; ++i)
	{
		for (int j = i; j < n; ++j)
		{
			Transposed.data[i][j] = data[j][i];
			Transposed.data[j][i] = data[i][j];
		}
		cout << '\n';
	}

	return Transposed;
}

void Matrix::Random()
{
	srand(time(0));
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			data[i][j] = 1 + rand() % 10;
}

Matrix Matrix::Pre_Minor(int row, int col) const //todo+
{
	if (n != m) throw Different_Dimensions();
	Matrix New_Matrix(m - 1, n - 1);
	int in = 0, jn = 0;

	for (int i = 0; i < m; i++)
	{
		if (i != row)
		{
			jn = 0;
			for (int j = 0; j < m; j++) {
				if (j != col) {
					New_Matrix.data[in][jn] = data[i][j];
					jn++;
				}
			}
			in++;
		}
	}

	return New_Matrix;
}

double Matrix::NDeterminant(int size)
{
	if (n != m) throw Different_Dimensions();
	Matrix TmpMatrix(m, m);
	int new_size = size - 1;
	double d = 0;
	int k = 1; //(-1) в степени i
	if (size < 1) cout << "Определитель вычислить невозможно!";
	if (m == 1)
	{
		d = data[0][0];
		return(d);
	}
	if (m == 2)
	{
		d = data[0][0] * data[1][1] - data[1][0] * data[0][1];
		return(d);
	}
	if (m > 2)
	{


		for (int i = 0; i < size; i++)
		{

			TmpMatrix = (*this).Pre_Minor(i, 0);

			d += k * data[i][0] * TmpMatrix.NDeterminant(new_size);
			k = -k;
		}

	}

	return(d);
}



Matrix Matrix::Search_Matrix_X(const Matrix& Vector)
{
	if (this->n != Vector.m) throw Different_Dimensions();

	double det = (*this).NDeterminant(m);
	if (det == 0) throw Zero_Determinant();

	cout << det << endl;
	Matrix Minors(n, m);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Minors.data[i][j] = (pow(-1, (i + j))) * Pre_Minor(i, j).NDeterminant(m - 1);
		}
	}
	cout << Minors << endl;

	Matrix Minors_Transpose = Minors.Transpose();
	cout << Minors_Transpose << endl;
	Matrix Ans = ((1 / det) * Minors_Transpose) * Vector;


	return Ans;
}

ostream& operator << (ostream& os, const Matrix& New_Matrix)
{
	for (int i = 0; i < New_Matrix.m; i++) {
		for (int j = 0; j < New_Matrix.n; j++) {
			os << "\t" << New_Matrix.Get_Data(i, j);
		}
		cout << endl;
	}
	return os;
}

