//
// Created by a.pesterev on 13.03.2023.
//

#ifndef POROUS_MEDIA_METHODS_HPP
#define POROUS_MEDIA_METHODS_HPP

constexpr FLOAT_TYPE two = FLOAT_TYPE(2);

__forceinline FLOAT_TYPE integrate(const std::function<FLOAT_TYPE(FLOAT_TYPE)>& function, const FLOAT_TYPE& a, const FLOAT_TYPE& b) {
    const FLOAT_TYPE dx = b - a;
    const FLOAT_TYPE avg_f = (function(b) + function(a)) / two;
    return dx * avg_f;
}

__forceinline FLOAT_TYPE interpolate(const FLOAT_TYPE& v_0, const FLOAT_TYPE& v_1, const FLOAT_TYPE& x_0, const FLOAT_TYPE& x_1, const FLOAT_TYPE& at_x) {
    const FLOAT_TYPE df = v_1 - v_0;
    const FLOAT_TYPE dx = x_1 - x_0;
    return (df/dx)*(at_x - x_0) + v_0;
}

__forceinline FLOAT_TYPE get_x(const INDEX_TYPE& i, const INDEX_TYPE& j, const std::vector<FLOAT_TYPE>& deltas) {
    FLOAT_TYPE x = 0;
    INTEGER_TYPE bix = i;
    INTEGER_TYPE biy = j;
    while (bix != 0) {
        x += deltas[biy];
        bix--;
        biy--;
    }
    return x;
}

#endif //POROUS_MEDIA_METHODS_HPP
