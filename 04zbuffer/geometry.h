#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// primary class template
template<size_t DimCols,size_t DimRows,typename T> class mat;

/**
 * @brief 向量类
 * 
 * @tparam DIM 向量的维度
 * @tparam T 向量的数据类型
 * 
 */
template <size_t DIM, typename T> struct vec
{
    vec()
    {
        // 初始化一个数组，数组的长度为DIM，数组中的每个元素都是T类型的默认值，0
        // 因为size_t是一个无符号的数，所以小于0的时候会变成一个很大的数，自动停止执行for循环
        for (size_t i=DIM; i--; data_[i] = T());
    }
    // 下标运算符重载，用于访问向量中的元素，一个用于非常量对象，一个用于常量对象
    T& operator[](const size_t i)
    {
        // assert是c++的一个宏，用于在运行时检查条件是否为真
        // 如果为假，输出错误信息，并终止程序。并返回失败的条件和代码的位置，方便调试
        assert(i<DIM);
        return data_[i];
    }
    
    const T& operator[](const size_t i)const
    {
        assert(i<DIM);
        return data_[i];
    }
private:
    T data_[DIM];
};

// 初始化二维向量
template <typename T> struct vec<2,T> {
    vec() : x(T()), y(T()) {}
    vec(T X, T Y) : x(X), y(Y) {}

    // 模板拷贝构造函数，可以从不同类型的vec对象拷贝数据
    template <class U> vec<2,T>(const vec<2,U> &v);
          T& operator[](const size_t i)       { assert(i<2); return i<=0 ? x : y; }
    const T& operator[](const size_t i) const { assert(i<2); return i<=0 ? x : y; }

    T x,y;
};

template <typename T> struct vec<3,T> {
    vec() : x(T()), y(T()), z(T()) {}
    vec(T X, T Y, T Z) : x(X), y(Y), z(Z) {}
    template <class U> vec<3,T>(const vec<3,U> &v);
          T& operator[](const size_t i)       { assert(i<3); return i<=0 ? x : (1==i ? y : z); }
    const T& operator[](const size_t i) const { assert(i<3); return i<=0 ? x : (1==i ? y : z); }
    
    // 模
    float norm() { return std::sqrt(x*x+y*y+z*z); }

    // 归一化
    // 函数的参数类型为T，默认值为1
    // this是一个指向当前对象的指针，*this则是当前对象的引用
    vec<3,T> & normalize(T l=1) { *this = (*this)*(l/norm()); return *this; }

    T x,y,z;
};

/////////////////////////////////////////////////////////////////////////////////
// 点乘
template<size_t DIM,typename T> T operator*(const vec<DIM,T>& lhs, const vec<DIM,T>& rhs) {
    T ret = T();        // 初始化默认为0
    for (size_t i=DIM; i--; ret+=lhs[i]*rhs[i]);
    return ret;
}
// 使用const修改参数：可以确保函数内部不会修改传入的参数
// 使用const引用可以避免拷贝临时对象，从而提高性能
// c++中，operator overloading是通过定义一个特殊的函数来实现的，而不是通过override关键字
// override通常用来重写基类的虚函数

// 如何判定是使用原始的*运算符还是重载的*运算符？
// 操作数的类型；模板匹配；作用域和可见性，确保运算符使用时是可见的


template<size_t DIM,typename T>vec<DIM,T> operator+(vec<DIM,T> lhs, const vec<DIM,T>& rhs) {
    for (size_t i=DIM; i--; lhs[i]+=rhs[i]);
    return lhs;
}

template<size_t DIM,typename T>vec<DIM,T> operator-(vec<DIM,T> lhs, const vec<DIM,T>& rhs) {
    for (size_t i=DIM; i--; lhs[i]-=rhs[i]);
    return lhs;
}

template<size_t DIM,typename T,typename U> vec<DIM,T> operator*(vec<DIM,T> lhs, const U& rhs) {
    for (size_t i=DIM; i--; lhs[i]*=rhs);
    return lhs;
}

template<size_t DIM,typename T,typename U> vec<DIM,T> operator/(vec<DIM,T> lhs, const U& rhs) {
    for (size_t i=DIM; i--; lhs[i]/=rhs);
    return lhs;
}

// 将低维向量嵌入到高维向量中
template<size_t LEN,size_t DIM,typename T> vec<LEN,T> embed(const vec<DIM,T> &v, T fill=1) {
    vec<LEN,T> ret;
    for (size_t i=LEN; i--; ret[i]=(i<DIM?v[i]:fill));
    return ret;
}

// 将高维DIM向量投影到低维LEN向量中
template<size_t LEN,size_t DIM, typename T> vec<LEN,T> proj(const vec<DIM,T> &v) {
    vec<LEN,T> ret;
    for (size_t i=LEN; i--; ret[i]=v[i]);
    return ret;
}

// TODO:叉乘推导
template <typename T> vec<3,T> cross(vec<3,T> v1, vec<3,T> v2) {
    return vec<3,T>(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

// 重载<<运算符，用于输出向量
template <size_t DIM, typename T> std::ostream& operator<<(std::ostream& out, vec<DIM,T>& v) {
    for(unsigned int i=0; i<DIM; i++) {
        out << v[i] << " " ;
    }
    return out ;
}

/////////////////////////////////////////////////////////////////////////////////

template<size_t DIM,typename T> struct dt {
    static T det(const mat<DIM,DIM,T>& src) {
        T ret=0;
        for (size_t i=DIM; i--; ret += src[0][i]*src.cofactor(0,i));
        return ret;
    }
};

template<typename T> struct dt<1,T> {
    static T det(const mat<1,1,T>& src) {
        return src[0][0];
    }
};

/////////////////////////////////////////////////////////////////////////////////

template<size_t DimRows,size_t DimCols,typename T> class mat {
    vec<DimCols,T> rows[DimRows];
public:
    mat() {}

    vec<DimCols,T>& operator[] (const size_t idx) {
        assert(idx<DimRows);
        return rows[idx];
    }

    const vec<DimCols,T>& operator[] (const size_t idx) const {
        assert(idx<DimRows);
        return rows[idx];
    }

    vec<DimRows,T> col(const size_t idx) const {
        assert(idx<DimCols);
        vec<DimRows,T> ret;
        for (size_t i=DimRows; i--; ret[i]=rows[i][idx]);
        return ret;
    }

    void set_col(size_t idx, vec<DimRows,T> v) {
        assert(idx<DimCols);
        for (size_t i=DimRows; i--; rows[i][idx]=v[i]);
    }

    static mat<DimRows,DimCols,T> identity() {
        mat<DimRows,DimCols,T> ret;
        for (size_t i=DimRows; i--; )
            for (size_t j=DimCols;j--; ret[i][j]=(i==j));
        return ret;
    }

    T det() const {
        return dt<DimCols,T>::det(*this);
    }

    mat<DimRows-1,DimCols-1,T> get_minor(size_t row, size_t col) const {
        mat<DimRows-1,DimCols-1,T> ret;
        for (size_t i=DimRows-1; i--; )
            for (size_t j=DimCols-1;j--; ret[i][j]=rows[i<row?i:i+1][j<col?j:j+1]);
        return ret;
    }

    T cofactor(size_t row, size_t col) const {
        return get_minor(row,col).det()*((row+col)%2 ? -1 : 1);
    }

    mat<DimRows,DimCols,T> adjugate() const {
        mat<DimRows,DimCols,T> ret;
        for (size_t i=DimRows; i--; )
            for (size_t j=DimCols; j--; ret[i][j]=cofactor(i,j));
        return ret;
    }

    mat<DimRows,DimCols,T> invert_transpose() {
        mat<DimRows,DimCols,T> ret = adjugate();
        T tmp = ret[0]*rows[0];
        return ret/tmp;
    }

    mat<DimRows,DimCols,T> invert() {
        return invert_transpose().transpose();
    }

    mat<DimCols,DimRows,T> transpose() {
        mat<DimCols,DimRows,T> ret;
        for (size_t i=DimCols; i--; ret[i]=this->col(i));
        return ret;
    }
};

/////////////////////////////////////////////////////////////////////////////////

template<size_t DimRows,size_t DimCols,typename T> vec<DimRows,T> operator*(const mat<DimRows,DimCols,T>& lhs, const vec<DimCols,T>& rhs) {
    vec<DimRows,T> ret;
    for (size_t i=DimRows; i--; ret[i]=lhs[i]*rhs);
    return ret;
}

template<size_t R1,size_t C1,size_t C2,typename T>mat<R1,C2,T> operator*(const mat<R1,C1,T>& lhs, const mat<C1,C2,T>& rhs) {
    mat<R1,C2,T> result;
    for (size_t i=R1; i--; )
        for (size_t j=C2; j--; result[i][j]=lhs[i]*rhs.col(j));
    return result;
}

template<size_t DimRows,size_t DimCols,typename T>mat<DimCols,DimRows,T> operator/(mat<DimRows,DimCols,T> lhs, const T& rhs) {
    for (size_t i=DimRows; i--; lhs[i]=lhs[i]/rhs);
    return lhs;
}

template <size_t DimRows,size_t DimCols,class T> std::ostream& operator<<(std::ostream& out, mat<DimRows,DimCols,T>& m) {
    for (size_t i=0; i<DimRows; i++) out << m[i] << std::endl;
    return out;
}

/////////////////////////////////////////////////////////////////////////////////

typedef vec<2,  float> Vec2f;
typedef vec<2,  int>   Vec2i;
typedef vec<3,  float> Vec3f;
typedef vec<3,  int>   Vec3i;
typedef vec<4,  float> Vec4f;
typedef mat<4,4,float> Matrix;
#endif //__GEOMETRY_H__

