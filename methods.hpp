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

#endif //POROUS_MEDIA_METHODS_HPP
