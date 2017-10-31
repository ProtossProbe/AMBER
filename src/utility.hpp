//
//  utility.hpp
//
//
//  Created by Protoss Probe on 2017/06/07.
//  Copyright © 2016-2017 probe. All rights reserved.
//

#ifndef _UTILITY_HPP_
#define _UTILITY_HPP_

#include <math.h>
namespace ProbeUtility {
template <class T, class S>
inline void joinVec(S &v, const T &v1, const T &v2) {
    int i = 0;
    for (auto var : v1) {
        v[i] = var;
        i++;
    }
    for (auto var : v2) {
        v[i] = var;
        i++;
    }
}

template <class T> inline double vecNorm(const T &v) {
    double result = 0;
    for (auto var : v) {
        result += var * var;
    }
    return sqrt(result);
}

template <class T> inline double vecDot(const T &v1, const T &v2) {
    double result = 0;
    for (int i = 0; i < v1.size(); i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

template <class T> inline void vec3Cross(T &v3, const T &v1, const T &v2) {
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
} // namespace ProbeUtility

#endif