#include "Matrix.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

Matrix::Matrix()
{
	rows = 0;
	cols = 0;
	data = NULL;
}

Matrix::Matrix(const Matrix& task) : rows(task.rows), cols(task.cols)
{
	try
	{
		data = new double[rows * cols];
	}
	catch (...)
	{
		cerr << "Memory allocation error" << endl;
		exit(1);
	}
	for (size_t i = 0; i < rows * cols; ++ i)
	{
		data[i] = task.data[i];
	}
}

Matrix::Matrix(size_t rows, size_t cols) : rows(rows), cols(cols)
{
	try
	{
		data = new double[rows * cols];
	}
	catch (...)
	{
		cerr << "Memory allocation error" << endl;
		exit(1);
	}
}

Matrix::~Matrix()
{
	delete[] data;
	data = NULL;
}

void Matrix::print() const
{
	for (size_t i = 0; i < rows; ++ i)
	{
		for (size_t j = 0; j < cols; ++ j)
		{
			cout << at(i, j) << " ";
		}
		cout << endl;
	}
}

void Matrix::print(const Matrix& task)
{
	for (size_t i = 0; i < task.rows; ++ i)
	{
		for (size_t j = 0; j < task.cols; ++ j)
		{
			cout << task.at(i, j) << " ";
		}
		cout << endl;
	}
}

double& Matrix::at (size_t r, size_t c)
{
	double& task = data[r * cols + c];
	return task;
}

const double& Matrix::at (size_t r, size_t c) const
{
	const double& task = data[r * cols + c];
	return task;
}

Matrix Matrix::add (const Matrix& l, const Matrix& r)
{
	Matrix task (l.rows, l.cols);
	
	for (size_t i = 0; i < task.rows; ++ i)
	{
		for (size_t j = 0; j < task.cols; ++ j)
			task.at(i, j) = l.at(i, j) + r.at(i, j);
	}
	return task;
}

Matrix Matrix::mul (const Matrix& l, const Matrix& r)
{
	
	Matrix task (l.rows, r.cols);
	
	for (size_t i = 0; i < task.rows; ++ i)
	{
		for (size_t j = 0; j < task.cols; ++ j)
		{
			task.at(i, j) = 0.0;
		}
	}
	for (size_t i = 0; i < task.rows; ++ i)
	{
		for (size_t j = 0; j < task.cols; ++ j)
		{
			for (size_t k = 0; k < l.cols; ++ k)
			{
				task.at(i, j) += l.at(i, k) * r.at(k, j);	
			}
		}
	}
	return task;
}

Matrix Matrix::mul (double a, const Matrix& m)
{
	
	Matrix task (m.rows, m.cols);
	
	for (size_t i = 0; i < task.rows; ++ i)
	{
		for (size_t j = 0; j < task.cols; ++ j)
		{
			task.at(i, j) = a * m.at(i, j);
		}
	}
	return task;
}

void Matrix::mul (double a)
{
	for (size_t i = 0; i < rows; ++ i)
	{
		for (size_t j = 0; j < cols; ++ j)
		{
			at(i, j) *= a;
		}
	}
}

void Matrix::Transpose ()
{
	size_t rows_result = cols;
	size_t cols_result = rows;
	double* data_result;
	
	try
	{
		data_result = new double[rows_result * cols_result];
	}
	catch (...)
	{
		cerr << "Memory allocation error" << endl;
		exit(1);
	}
	
	for (size_t i = 0; i < rows_result; ++ i)
	{
		for (size_t j = 0; j < cols_result; ++ j)
		{
			data_result[i * cols_result + j] = at(j, i);
		}
	}
	
	delete[] data;
	
	data = data_result;
	rows = rows_result;
	cols = cols_result;
}

Matrix Matrix::operator +(const Matrix& b) const
{
    return add(*this, b);
}

Matrix Matrix::operator -(const Matrix& b) const
{
    double x = -1;
    return add(*this, mul(x,b));
}

Matrix Matrix::operator *(const Matrix& b) const
{
    return mul(*this, b);
}

Matrix Matrix::operator *(double d) const
{
    return mul(d, *this);
}

Matrix operator *(double d, Matrix& m)
{
    return Matrix::mul(d, m);
}

std::istream& operator >>(std::istream& input_m, Matrix& matrix)
{
    input_m >> matrix.rows >> matrix.cols;
    try
	{
		matrix.data = new double[matrix.rows*matrix.cols];
	}
	catch (...)
	{
		cerr << "Erorr allocating memory" << endl;
		exit(1);
	}
    for (size_t i = 0; i < matrix.rows; ++i)
    {
		for (size_t j = 0; j < matrix.cols; ++j)
		{
		    double cur;
		    input_m >> cur;
		    matrix.at(i, j) = cur;
		}
	}
    return input_m;
}

std::ostream& operator <<(std::ostream& output_m, const Matrix& matrix)
{
    output_m << matrix.rows << " " << matrix.cols << endl;
    for (size_t i = 0; i < matrix.rows; ++i)
    {
		for (size_t j = 0; j < matrix.cols; ++j)
		{
		    output_m << matrix.at(i, j) << " ";
		}
		output_m << endl;
	}
	return output_m;
}

Matrix Matrix::operator =(const Matrix& matrix)
{
    this->rows = matrix.rows;
    this->cols = matrix.cols;

	try
	{
		this->data = new double[this->rows*this->cols];
	}
	catch (...)
	{
		cerr << "Erorr allocating memory" << endl;
		exit(1);
	}
    for (size_t i = 0; i < this->rows; ++i)
        for (size_t j = 0; j < this->cols; ++j)
        {
            this->at(i,j) = matrix.at(i,j);
        }
    return *this;
}

Matrix Matrix::Upper() const
{
    Matrix task;
    task = *this;
    for (size_t i = 0; i < std::min(task.rows, task.cols); ++i)
    {
        size_t task_max = i;
        double another_task_max = task.at(i,i);
        for (size_t j = i+1; j < std::min(task.rows, task.cols); ++j)
        {
            if (fabs(task.at(j,i)) > fabs(another_task_max))
            {
                task_max = j;
                another_task_max = task.at(j,i);
            }
        }
        if (fabs(another_task_max) <= 0.0001)
            continue;
        for (size_t j = i; j < task.cols; ++j)
            std::swap (task.at(i, j), task.at(task_max, j));
        for (size_t j = i+1; j < std::min(task.rows, task.cols); ++j)
        {
            for (size_t x = i+1; x < task.cols; ++x)
            {
                task.at(j, x) -= task.at(i, x)*task.at(j,i)/task.at(i,i);
            }
            task.at(j, i) -= task.at(i, i)*task.at(j,i)/task.at(i,i);
        }
    }
    return task;
}

bool Matrix::is_triangular() const
{

    bool task = true;
    for (size_t i = 0; i < rows && task; ++i)
    {
        for (size_t j = 0; j < i; ++j)
            if (at(i,j) != 0)
            {
                task = false;
                break;
            }
    }
    if (task)
        return task;
    task = true;
    for (size_t i = 0; i < rows && task; ++i)
    {
        for (size_t j = i+1; j < cols; ++j)
            if (at(i,j) != 0)
            {
                task = false;
                break;
            }
    }
    return task;
}