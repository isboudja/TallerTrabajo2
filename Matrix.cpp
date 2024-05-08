#include "Matrix.h"
#include <iostream>
#include <iomanip>
#include "math.h"

Matrix::Matrix(int fil, int col)
{
    this->fil = fil;
    this->col = col;
    initMatrix();
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
                matrix[i][j] = 0;
        }
}
 
Matrix::Matrix(int fil, int col, double v[], int n)
{


    this->fil = fil;
    this->col = col;
    initMatrix();


    int k = 0;
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}
 
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}

int Matrix::fils()
{
    return fil;
}

int Matrix::cols()
{
    return col;
}

Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}
 
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}
 
Matrix& Matrix::operator=(const Matrix& matrix2)
{
        if (this != &matrix2) {
            this->fil = matrix2.fil;
            this->col = matrix2.col;
            for (int i = 0; i < fil; i++) {
                for (int j = 0; j < col; j++) {
                    this->matrix[i][j] = matrix2.matrix[i][j];
                }
            }
        }
        return *this;
    }

Matrix Matrix::operator^(double exponent) const {
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            result.matrix[i][j] = pow(matrix[i][j], exponent);
        }
    }

    return result;
}
 
Matrix Matrix::operator+(const Matrix& matrix2) const
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator-(const Matrix& matrix2) const
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}

Matrix Matrix::operator/(const Matrix& matrix2) const {
    if (fil != matrix2.fil || col != matrix2.col) {
        return Matrix(fil, col);
    }

    Matrix result(fil, col);

    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            if (matrix2.matrix[i][j] != 0) {
                result.matrix[i][j] = matrix[i][j] / matrix2(i + 1, j + 1);
            } else {
                result.matrix[i][j] = 0;
            }
        }
    }

    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            result.matrix[i][j] = matrix[i][j] * scalar;
        }
    }

    return result;
}

double Matrix::norm() const {
    double sumOfSquares = 0.0;

    if (col < 3)
        throw "Error";

    for (int i = 0; i < this->col; ++i)
        sumOfSquares += matrix[0][i] * matrix[0][i];

    return sqrt(sumOfSquares);
}

Matrix Matrix::operator/(double scalar) const {
    Matrix result(fil, col);
    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            result(i + 1, j + 1) = matrix[i][j] / scalar;
        }
    }
    return result;
}
 
Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(fil, col);
 
    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }
 
    return result;
}
 
 
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}
 
void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}