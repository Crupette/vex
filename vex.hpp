#ifndef VEX_H
#define VEX_H 1

#include <array>
#include <cmath>
#include <cstddef>
#include <ostream>

namespace vex
{

template<typename T, unsigned D>
struct vec_dimd
{
    std::array<T, D> v;
    
    constexpr vec_dimd() = default;
    template<typename ...Args>
    constexpr vec_dimd(Args&& ...args) : v{args...} {}
    constexpr vec_dimd(T args[D]) : v(args) {}
    constexpr vec_dimd(T fill) { v.fill(fill); }
    constexpr vec_dimd(const vec_dimd<T, D-1> &vec, T val) { std::copy(vec.v.begin(), vec.v.end(), v); v[D-1] = val; }

    constexpr T operator[](std::size_t i) { return v[i]; }

    vec_dimd<T, D> operator+=(const vec_dimd<T, D> &rhs) { for(size_t i = 0; i < D; i++) v[i] += rhs.v[i]; return *this; }
    vec_dimd<T, D> operator-=(const vec_dimd<T, D> &rhs) { for(size_t i = 0; i < D; i++) v[i] -= rhs.v[i]; return *this; }
    vec_dimd<T, D> operator*=(const vec_dimd<T, D> &rhs) { for(size_t i = 0; i < D; i++) v[i] *= rhs.v[i]; return *this; }
    vec_dimd<T, D> operator/=(const vec_dimd<T, D> &rhs) { for(size_t i = 0; i < D; i++) v[i] /= rhs.v[i]; return *this; }

    vec_dimd<T, D> operator*=(const T &rhs) { for(size_t i = 0; i < D; i++) v[i] *= rhs; return *this; }
    vec_dimd<T, D> operator/=(const T &rhs) { for(size_t i = 0; i < D; i++) v[i] *= rhs; return *this; }

    auto operator<=>(const vec_dimd<T,D> &rhs) const { return magnitude() <=> rhs.magnitude(); }
    bool operator==(const vec_dimd<T,D> &rhs) const = default;

    friend vec_dimd<T, D> operator+(vec_dimd<T, D> lhs, const vec_dimd<T, D> &rhs) { lhs += rhs; return lhs; }
    friend vec_dimd<T, D> operator-(vec_dimd<T, D> lhs, const vec_dimd<T, D> &rhs) { lhs -= rhs; return lhs; }
    friend vec_dimd<T, D> operator*(vec_dimd<T, D> lhs, const vec_dimd<T, D> &rhs) { lhs *= rhs; return lhs; }
    friend vec_dimd<T, D> operator/(vec_dimd<T, D> lhs, const vec_dimd<T, D> &rhs) { lhs /= rhs; return lhs; }

    friend vec_dimd<T, D> operator*(vec_dimd<T, D> lhs, const T &rhs) { lhs *= rhs; return lhs; }
    friend vec_dimd<T, D> operator/(vec_dimd<T, D> lhs, const T &rhs) { lhs /= rhs; return lhs; }    

    friend std::ostream &operator<<(std::ostream &os, const vec_dimd<T,D> obj) {
        os << '{';
        for(std::size_t i = 0; i < D; i++) os << obj.v[i] << ((i < (D - 1)) ? ',' : '}');
        return os;
    }

    T magnitude() const { T t; for(size_t i = 0; i < D; i++) t += v[i] * v[i]; return sqrt(t); }
    T dot(const vec_dimd<T, D> &b) const { T t = v[0] * b.v[0]; for(size_t i = 1; i < D; i++) t += (v[i] * b.v[i]); return t; }
};

template<typename T>
using vec2 = vec_dimd<T, 2>;

template<typename T>
using vec3 = vec_dimd<T, 3>;

};

#endif
