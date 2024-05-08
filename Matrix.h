#ifndef _MATRIX_
#define _MATRIX_

class Matrix
{
    public:
        Matrix(int fil, int col);
        Matrix(int fil, int col, double v[], int n);
        int fils();
        int cols();
        int fil;
        int col;
        Matrix(const Matrix& m);
        ~Matrix();
 
        Matrix& operator=(const Matrix& matrix2);
        Matrix  operator+(const Matrix& matrix2) const;
        Matrix  operator-(const Matrix& matrix2) const;
        Matrix operator^(double exponent) const;
        Matrix  operator*(const Matrix& matrix2);
        Matrix operator/(const Matrix& matrix2) const;
        Matrix operator*(double scalar) const;
        Matrix operator/(double scalar) const;
    double& operator()(const int i, const int j) const;
        double norm() const;


        void print();
 
    private:
        void initMatrix();
 
    private:

        double **matrix;
};

#endif
