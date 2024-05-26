#ifndef _MATRIX_
#define _MATRIX_

class Matrix
{
    public:
        Matrix(int fil, int col);
        Matrix(int fil, int col, double v[], int n);

        int fils();
        int cols();
        Matrix sub(int i,int j) const;
        static Matrix concat(const Matrix& mat1, const Matrix& mat2);
        int fil;
        int col;
        Matrix(const Matrix& m);
        ~Matrix();
        Matrix transpose() const;
        Matrix& operator=(const Matrix& matrix2);
        Matrix operator+(double scalar) const;
        Matrix  operator+(const Matrix& matrix2) const;
        Matrix  operator-(const Matrix& matrix2) const;
        Matrix operator^(double exponent) const;
        Matrix  operator*(const Matrix& matrix2);
        Matrix operator/(const Matrix& matrix2) const;
        Matrix operator*(double scalar) const;
        Matrix operator/(double scalar) const;
        double& operator()(const int i, const int j) const;
        double norm() const;
        Matrix operator-(double scalar) const;
        void print();
 
    private:
        void initMatrix();
 
    private:

        double **matrix;
};

#endif
