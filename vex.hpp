#ifndef VEX_H
#define VEX_H 1

#include <array>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <concepts>

namespace vex
{

template<typename T>
concept vectype = std::is_integral<T>::value || std::is_floating_point<T>::value;

template<vectype T, unsigned D>
    requires (D > 1)
struct vec_dimd
{
    std::array<T, D> v;
    
    explicit vec_dimd() = default;
    template<typename ...Args>
    explicit vec_dimd(Args&&... args) : v{args...} {}
    explicit vec_dimd(T args[D]) : v(args) {}
    explicit vec_dimd(T fill) { v.fill(fill); }

    T operator[](std::size_t i) { return v[i]; }
    T operator[](std::size_t i) const { return v[i]; }

    vec_dimd<T, D> operator+=(const vec_dimd<T, D> &rhs) { for(size_t i = 0; i < D; i++) v[i] += rhs.v[i]; return *this; }
    vec_dimd<T, D> operator-=(const vec_dimd<T, D> &rhs) { for(size_t i = 0; i < D; i++) v[i] -= rhs.v[i]; return *this; }
    vec_dimd<T, D> operator*=(const vec_dimd<T, D> &rhs) { for(size_t i = 0; i < D; i++) v[i] *= rhs.v[i]; return *this; }
    vec_dimd<T, D> operator/=(const vec_dimd<T, D> &rhs) { for(size_t i = 0; i < D; i++) v[i] /= rhs.v[i]; return *this; }

    vec_dimd<T, D> operator*=(const T &rhs) { for(size_t i = 0; i < D; i++) v[i] *= rhs; return *this; }
    vec_dimd<T, D> operator/=(const T &rhs) { for(size_t i = 0; i < D; i++) v[i] /= rhs; return *this; }

    auto operator<=>(const vec_dimd<T,D> &rhs) const { return magnitude() <=> rhs.magnitude(); }
    bool operator==(const vec_dimd<T,D> &rhs) const { for(unsigned i = 0; i < D; i++) if(v[i] != rhs.v[i]) return false; return true; }

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

    /*Finds the distance from the origin to the ray cast in D dimension space using components of vec*/
    T magnitude() const { T t{}; for(size_t i = 0; i < D; i++) t += v[i] * v[i]; return std::sqrt(t); }

    /*Finds the dot product of itself and another vec of same T and D*/
    T dot(const vec_dimd<T, D> &b) const { T t; for(size_t i = 0; i < D; i++) t += (v[i] * b.v[i]); return t; }
    friend T dot(const vec_dimd<T,D> &lhs, const vec_dimd<T,D> &rhs) { return lhs.dot(rhs); }

    auto cross(const vec_dimd<T,D> &o) const {
        if constexpr(D == 2) return v[0]*o[0] - v[1]*o[1];
        if constexpr(D == 3) return vec_dimd<T, D>{v[1]*o[2]-v[2]*o[1], v[2]*o[0]-v[0]*o[2],v[0]*o[1]-v[1]*o[0]};
    }
};

template<vectype T>
using vec2 = vec_dimd<T, 2>;

template<vectype T>
using vec3 = vec_dimd<T, 3>;

/*(x, y) -> (r, theta)*/
template<vectype R, vectype P>
static vec2<R> 
polar(const vec2<P> &in)
{
    return vec2<R>(
                std::sqrt(in[0] * in[0] + in[1] * in[1]),
                std::atan2(in[1], in[0])
            );
}
/*(r, theta) -> (x, y)*/
template<vectype R, vectype P>
static vec2<R> 
cartesian(const vec2<P> &in)
{
    return vec2<R>(
                cos(in[1]) * in[0],
                sin(in[1]) * in[0]
            );
}

};

#endif
