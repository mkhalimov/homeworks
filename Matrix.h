#ifndef Matrix_h
#define Matrix_h 20221115L
#define MATRIX_HAS_GAUSS

#include <cstddef>
#include <istream>
#include <ostream>

class Matrix {
	 size_t rows;
	 size_t cols;
	 
	 double *data;
	 
	public:
	 
	 Matrix();
	 Matrix(const Matrix&);
	 Matrix(size_t rows, size_t cols);
	~Matrix();
	 
	 void print() const;
	 static void print(const Matrix&);
	 double& at(size_t r, size_t c);
	 const double& at(size_t r, size_t c) const ;
	 
	 static Matrix add(const Matrix& l, const Matrix& r);
	 static Matrix mul(const Matrix& l, const Matrix& r);
	 static Matrix mul(double a, const Matrix& m);
	 
	 void mul(double a);
	 void Transpose();
	 Matrix operator =(const Matrix&);
	 Matrix operator +(const Matrix&) const;
	 Matrix operator -(const Matrix&) const;
	 Matrix operator *(const Matrix&) const;
	 Matrix operator *(double) const;
     
	 friend Matrix operator *(double, const Matrix&);
	 friend std::istream & operator >>(std::istream&, Matrix&);
	 friend std::ostream& operator <<(std::ostream& out, const Matrix&);
     #ifdef MATRIX_HAS_GAUSS
        Matrix Upper() const;
        bool is_triangular() const;
     #endif /*MATRIX_HAS_GAUSS*/
};

#endif /*Matrix_h*/
