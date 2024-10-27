#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include "geometry.h"

// TODO:特化构造函数
// template <> template <> 是模板特化的语法，表示正在定义一个特化版本的模板函数
// 定义了两个模板类 Vec3 的特化构造函数，用于在不同的类型之间转接。
// 将 Vec3<float> 类型的对象转换为 Vec3<int> 类型的对象。
template <> template <> Vec3<int>::Vec3(const Vec3<float> &v) : x(int(v.x+.5)), y(int(v.y+.5)), z(int(v.z+.5)) {
    // 使用+0.5，然后强制转换为int，实现四舍五入
}

template <> template <> Vec3<float>::Vec3(const Vec3<int> &v) : x(v.x), y(v.y), z(v.z) {
}

// 注意在旧版本的 C++ 中（C++11 之前），在两个右尖括号之间需要有一个空格，以避免与 >> 操作符混淆。
// 但是在 C++11 中，这种情况已经被修复了，所以在 C++11 中，两个右尖括号之间不需要有空格。
// 初始化一个二维向量，外层向量被初始化为r个元素，每个元素都是一个浮点数向量。每个内层向量又被初始化为c个元素，每个元素都是0
Matrix::Matrix(int r, int c) : m(std::vector<std::vector<float> >(r, std::vector<float>(c, 0.f))), rows(r), cols(c) { }
// 一共初始化三个成员变量，m是一个二维向量，rows和cols分别是矩阵的行数和列数
// m被初始化为一个r行c列的二维向量，每个元素都是0

int Matrix::nrows() {
    return rows;
}

int Matrix::ncols() {
    return cols;
}

Matrix Matrix::identity(int dimensions) {
    Matrix E(dimensions, dimensions);
    for (int i=0; i<dimensions; i++) {
        for (int j=0; j<dimensions; j++) {
            E[i][j] = (i==j ? 1.f : 0.f);
        }
    }
    return E;
}

std::vector<float>& Matrix::operator[](const int i) {
    assert(i>=0 && i<rows);
    return m[i];
}

// NOTE:cols是函数的成员变量，因此在成员函数中可以直接使用，而无需声明
Matrix Matrix::operator*(const Matrix& a) {
    assert(cols == a.rows);
    Matrix result(rows, a.cols);
    for (int i=0; i<rows; i++) {
        for (int j=0; j<a.cols; j++) {
            result.m[i][j] = 0.f;
            for (int k=0; k<cols; k++) {
                result.m[i][j] += m[i][k]*a.m[k][j];
            }
        }
    }
    return result;
}

// 转置矩阵
Matrix Matrix::transpose() {
    Matrix result(cols, rows);
    for(int i=0; i<rows; i++)
        for(int j=0; j<cols; j++)
            result[j][i] = m[i][j];
    return result;
}


// 逆矩阵
Matrix Matrix::inverse() {
    assert(rows==cols);
    // augmenting the square matrix with the identity matrix of the same dimensions a => [ai]
    Matrix result(rows, cols*2);
    for(int i=0; i<rows; i++)
        for(int j=0; j<cols; j++)
            result[i][j] = m[i][j];
    for(int i=0; i<rows; i++)
        result[i][i+cols] = 1;
    // first pass
    for (int i=0; i<rows-1; i++) {
        // normalize the first row
        for(int j=result.cols-1; j>=0; j--)
            result[i][j] /= result[i][i];
        for (int k=i+1; k<rows; k++) {
            float coeff = result[k][i];
            for (int j=0; j<result.cols; j++) {
                result[k][j] -= result[i][j]*coeff;
            }
        }
    }
    // normalize the last row
    for(int j=result.cols-1; j>=rows-1; j--)
        result[rows-1][j] /= result[rows-1][rows-1];
    // second pass
    for (int i=rows-1; i>0; i--) {
        for (int k=i-1; k>=0; k--) {
            float coeff = result[k][i];
            for (int j=0; j<result.cols; j++) {
                result[k][j] -= result[i][j]*coeff;
            }
        }
    }
    // cut the identity matrix back
    Matrix truncate(rows, cols);
    for(int i=0; i<rows; i++)
        for(int j=0; j<cols; j++)
            truncate[i][j] = result[i][j+cols];
    return truncate;
}

std::ostream& operator<<(std::ostream& s, Matrix& m) {
    for (int i=0; i<m.nrows(); i++)  {
        for (int j=0; j<m.ncols(); j++) {
            s << m[i][j];
            if (j<m.ncols()-1) s << "\t";
        }
        s << "\n";
    }
    return s;
}

