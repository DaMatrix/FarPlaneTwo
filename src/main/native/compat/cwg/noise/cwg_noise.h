/*
 * Adapted from The MIT License (MIT)
 *
 * Copyright (c) 2020-2021 DaPorkchop_
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
 * is furnished to do so, subject to the following conditions:
 *
 * Any persons and/or organizations using this software must include the above copyright notice and this permission notice,
 * provide sufficient credit to the original authors of the project (IE: DaPorkchop_), as well as provide a link to the original project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#ifndef CWG_NOISE_H
#define CWG_NOISE_H

#include <fp2.h>
#include <fp2/simd.h>

#include <cassert>
#include <immintrin.h>

namespace fp2::cwg::noise {
    /**
     * An array of random vectors used for noise generation.
     *
     * These are stored as floats rather than doubles in order to minimize their impact on the CPU cache.
     */
    inline float RANDOM_VECTORS[1024];

    inline void setRandomVectors(JNIEnv* env, jfloatArray in) {
        int32_t in_length = env->GetArrayLength(in);
        if (in_length * sizeof(float) != sizeof(RANDOM_VECTORS)) {
            throw fp2::error("invalid array length", in_length);
        }

        env->GetFloatArrayRegion(in, 0, in_length, &RANDOM_VECTORS[0]);
    }

    constexpr double raw_to_double(uint64_t in) {
        const union {
            uint64_t i;
            double d;
        } u = {in};
        return u.d;
    }

    template<size_t VEC_LANES> inline typename fp2::simd::type_vec<double, VEC_LANES>::TYPE makeInt32Range(
            typename fp2::simd::type_vec<double, VEC_LANES>::TYPE n) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;

        DOUBLE nn = n + n;
        DOUBLE magic0 = raw_to_double(0x3E10000000000000);
        DOUBLE magic1 = raw_to_double(0x41D0000000000000);

        return select(abs(n) >= 1073741824.0, nn - truncate(nn * magic0) * magic1 - sign_combine(magic1, n), n);

        // equivalent code:
        /*double d0 = n, d1, d2;
        if (n >= (d1 = 1073741824.0)) {
            d0 = d0 + d0;
            d1 = raw_to_double(0x3E10000000000000);
            d1 = d0 * d1;
            d1 = trunc(d1); // frintz d1, d1
            d2 = raw_to_double(0x41D0000000000000);
            d1 = d1 * d2;
            d0 = d0 - d1;
            d0 = d0 - d2;
        } else if (n <= (d1 = -1073741824.0)) {
            d0 = d0 + d0;
            d1 = raw_to_double(0x3E10000000000000);
            d1 = d0 * d1;
            d1 = trunc(d1); //frintz d1, d1
            d2 = raw_to_double(0x41D0000000000000);
            d1 = d1 * d2;
            d0 = d0 - d1;
            d0 = d0 + d2;
        }
        return d0;*/

        //equivalent code:
        /*if (n >= 1073741824.0) {
            return (fmod(2.0 * n, 1073741824.0)) - 1073741824.0;
        } else if (n <= -1073741824.0) {
            return (fmod(2.0 * n, 1073741824.0)) + 1073741824.0;
        } else {
            return n;
        }*/
    }

    template<typename T> constexpr T min(T a, T b) {
        return a < b ? a : b;
    }

    template<typename T> constexpr T max(T a, T b) {
        return a > b ? a : b;
    }

    template<typename T> constexpr T clamp(T val, T min, T max) {
        return val > min ? val < max ? val : max : min;
    }

    template<typename T> constexpr T abs(T a) {
        return a < 0 ? -a : a;
    }

    template<typename T> constexpr T copysign(T src, T dst) {
        return src < 0 ? -abs(dst) : abs(dst);
    }

    template<typename T> constexpr T signum(T n) {
        return n < 0 ? -1 : n > 0 ? 1 : 0;
    }

    template<typename T> constexpr T sCurve3(T a) {
        return a * a * (3.0 - 2.0 * a);
    }

    template<typename T> constexpr T vectorIndex(T ix, T iy, T iz, T seed) {
        T vectorIndex = 1619 * ix + 31337 * iy + 6971 * iz + 1013 * seed;
        return (vectorIndex ^ (vectorIndex >> 8)) & 0xFF;
    }

    template<typename T> constexpr T lerp(T n0, T n1, T a) {
        return n0 + (n1 - n0) * a;
    }

    template<size_t VEC_LANES> inline typename fp2::simd::type_vec<double, VEC_LANES>::TYPE gradientNoise3d(
            typename fp2::simd::type_vec<double, VEC_LANES>::TYPE fx, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE fy, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE fz,
            typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE ix, typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE iy, typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE iz,
            typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE seed) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;
        using INT = typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE;

        INT vi = vectorIndex(ix, iy, iz, seed) << 2;
        DOUBLE xvGradient, yvGradient, zvGradient;

        if constexpr (VEC_LANES == 8 && fp2::simd::LANES_32 == 16 && fp2::simd::LANES_64 == 8) { //AVX512
            using FLOAT4 = typename fp2::simd::type_vec<double, 4>::FLOAT::TYPE;
            using FLOAT8 = typename fp2::simd::type_vec<double, 8>::FLOAT::TYPE;
            using FLOAT16 = typename fp2::simd::type_vec<float, 16>::TYPE;

            //load 8 SSE registers in XYZW XYZW XYZW XYZW XYZW XYZW XYZW XYZW format
            FLOAT4 f0 = FLOAT4().load(&RANDOM_VECTORS[vi[0]]);
            FLOAT4 f1 = FLOAT4().load(&RANDOM_VECTORS[vi[1]]);
            FLOAT4 f2 = FLOAT4().load(&RANDOM_VECTORS[vi[2]]);
            FLOAT4 f3 = FLOAT4().load(&RANDOM_VECTORS[vi[3]]);
            FLOAT4 f4 = FLOAT4().load(&RANDOM_VECTORS[vi[4]]);
            FLOAT4 f5 = FLOAT4().load(&RANDOM_VECTORS[vi[5]]);
            FLOAT4 f6 = FLOAT4().load(&RANDOM_VECTORS[vi[6]]);
            FLOAT4 f7 = FLOAT4().load(&RANDOM_VECTORS[vi[7]]);

            //concatenate to get 4 AVX registers in XYZWXYZW XYZWXYZW XYZWXYZW XYZWXYZW format
            FLOAT8 f01 = FLOAT8(f0, f1);
            FLOAT8 f23 = FLOAT8(f2, f3);
            FLOAT8 f45 = FLOAT8(f4, f5);
            FLOAT8 f67 = FLOAT8(f6, f7);

            //concatenate to get 2 AVX512 registers in XYZWXYZWXYZWXYZW XYZWXYZWXYZWXYZW format
            FLOAT16 f0123 = FLOAT16(f01, f23);
            FLOAT16 f4567 = FLOAT16(f45, f67);

            //shuffle to get 3 AVX512 registers in XXXXXXXX???????? YYYYYYYY???????? ZZZZZZZZ???????? format
            FLOAT16 fx01234567________ = blend16<0, 4, 8, 12, 16, 20, 24, 28, V_DC, V_DC, V_DC, V_DC, V_DC, V_DC, V_DC, V_DC>(f0123, f4567);
            FLOAT16 fy01234567________ = blend16<1, 5, 9, 13, 17, 21, 25, 29, V_DC, V_DC, V_DC, V_DC, V_DC, V_DC, V_DC, V_DC>(f0123, f4567);
            FLOAT16 fz01234567________ = blend16<2, 6, 10, 14, 18, 22, 26, 30, V_DC, V_DC, V_DC, V_DC, V_DC, V_DC, V_DC, V_DC>(f0123, f4567);

            //truncate to get 3 AVX registers in XXXXXXXX YYYYYYYY ZZZZZZZZ format, then extend to double to get 3 AVX512 registers
            xvGradient = to_double(fx01234567________.get_low());
            yvGradient = to_double(fy01234567________.get_low());
            zvGradient = to_double(fz01234567________.get_low());
        } else if constexpr (VEC_LANES == 4 && fp2::simd::LANES_32 == 8 && fp2::simd::LANES_64 == 4) { //AVX and AVX2
            using FLOAT4 = typename fp2::simd::type_vec<double, 4>::FLOAT::TYPE;
            using FLOAT8 = typename fp2::simd::type_vec<double, 8>::FLOAT::TYPE;

            //load 4 SSE registers in XYZW XYZW XYZW XYZW format
            FLOAT4 f0 = FLOAT4().load(&RANDOM_VECTORS[vi[0]]);
            FLOAT4 f1 = FLOAT4().load(&RANDOM_VECTORS[vi[1]]);
            FLOAT4 f2 = FLOAT4().load(&RANDOM_VECTORS[vi[2]]);
            FLOAT4 f3 = FLOAT4().load(&RANDOM_VECTORS[vi[3]]);

            //concatenate to get 2 AVX registers in XYZWXYZW XYZWXYZW format
            FLOAT8 f01 = FLOAT8(f0, f1);
            FLOAT8 f23 = FLOAT8(f2, f3);

            //shuffle to get 3 AVX registers in XXXX???? YYYY???? ZZZZ???? format
            FLOAT8 fx0123____ = blend8<0, 4, 8, 12, V_DC, V_DC, V_DC, V_DC>(f01, f23);
            FLOAT8 fy0123____ = blend8<1, 5, 9, 13, V_DC, V_DC, V_DC, V_DC>(f01, f23);
            FLOAT8 fz0123____ = blend8<2, 6, 10, 14, V_DC, V_DC, V_DC, V_DC>(f01, f23);

            //truncate to get 3 SSE registers in XXXX YYYY ZZZZ format, then extend to double to get 3 AVX registers
            xvGradient = to_double(fx0123____.get_low());
            yvGradient = to_double(fy0123____.get_low());
            zvGradient = to_double(fz0123____.get_low());
        } else if constexpr (VEC_LANES == 4 && fp2::simd::LANES_32 == 4 && fp2::simd::LANES_64 == 2) { //SSE (with emulated 256-bit vectors)
            using FLOAT4 = typename fp2::simd::type_vec<double, 4>::FLOAT::TYPE;

            //load 4 SSE registers in XYZW XYZW XYZW XYZW format
            FLOAT4 f0 = FLOAT4().load(&RANDOM_VECTORS[vi[0]]);
            FLOAT4 f1 = FLOAT4().load(&RANDOM_VECTORS[vi[1]]);
            FLOAT4 f2 = FLOAT4().load(&RANDOM_VECTORS[vi[2]]);
            FLOAT4 f3 = FLOAT4().load(&RANDOM_VECTORS[vi[3]]);

            //shuffle to get 4 SSE registers in XXYY ZZWW XXYY ZZWW format
            FLOAT4 f01_XY = blend4<0, 4, 1, 5>(f0, f1);
            FLOAT4 f01_ZW = blend4<2, 6, V_DC, V_DC>(f0, f1);
            FLOAT4 f23_XY = blend4<0, 4, 1, 5>(f2, f3);
            FLOAT4 f23_ZW = blend4<2, 6, V_DC, V_DC>(f2, f3);

            //shuffle to get 3 SSE registers in XXXX YYYY ZZZZ format, then extend to double (emulated 256-bit, 2x SSE registers per variable)
            xvGradient = to_double(blend4<0, 1, 4, 5>(f01_XY, f23_XY));
            yvGradient = to_double(blend4<2, 3, 6, 7>(f01_XY, f23_XY));
            zvGradient = to_double(blend4<0, 1, 4, 5>(f01_ZW, f23_ZW));
        } else {
            assert(false);
        }

        return xvGradient * fx + yvGradient * fy + zvGradient * fz + 0.5;
    }

    template<size_t VEC_LANES> inline typename fp2::simd::type_vec<double, VEC_LANES>::TYPE noise3d(
            typename fp2::simd::type_vec<double, VEC_LANES>::TYPE x, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE y, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE z,
            typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE seed) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;
        using INT = typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE;

        //floor coordinates
        DOUBLE floorX = if_sub(x == 0.0, floor(x), 1.0);
        DOUBLE floorY = if_sub(y == 0.0, floor(y), 1.0);
        DOUBLE floorZ = if_sub(z == 0.0, floor(z), 1.0);

        //convert floored coordinates to ints to get the integer coordinates at the corners
        INT x0 = round_to_int32(floorX);
        INT x1 = x0 + 1;
        INT y0 = round_to_int32(floorY);
        INT y1 = y0 + 1;
        INT z0 = round_to_int32(floorZ);
        INT z1 = z0 + 1;

        //get coordinate's fractional parts
        DOUBLE fx = x - floorX;
        DOUBLE fy = y - floorY;
        DOUBLE fz = z - floorZ;

        //compute gradient vectors at each corner
        DOUBLE n000 = gradientNoise3d<VEC_LANES>(fx, fy, fz, x0, y0, z0, seed);
        DOUBLE n001 = gradientNoise3d<VEC_LANES>(fx, fy, fz - 1.0, x0, y0, z1, seed);
        DOUBLE n010 = gradientNoise3d<VEC_LANES>(fx, fy - 1.0, fz, x0, y1, z0, seed);
        DOUBLE n011 = gradientNoise3d<VEC_LANES>(fx, fy - 1.0, fz - 1.0, x0, y1, z1, seed);
        DOUBLE n100 = gradientNoise3d<VEC_LANES>(fx - 1.0, fy, fz, x1, y0, z0, seed);
        DOUBLE n101 = gradientNoise3d<VEC_LANES>(fx - 1.0, fy, fz - 1.0, x1, y0, z1, seed);
        DOUBLE n110 = gradientNoise3d<VEC_LANES>(fx - 1.0, fy - 1.0, fz, x1, y1, z0, seed);
        DOUBLE n111 = gradientNoise3d<VEC_LANES>(fx - 1.0, fy - 1.0, fz - 1.0, x1, y1, z1, seed);

        //smooth fractional coordinates
        DOUBLE xs = sCurve3(fx);
        DOUBLE ys = sCurve3(fy);
        DOUBLE zs = sCurve3(fz);

        //lerp between corners to get a single output sample
        DOUBLE n000_n100 = lerp(n000, n100, xs);
        DOUBLE n001_n101 = lerp(n001, n101, xs);
        DOUBLE n010_n110 = lerp(n010, n110, xs);
        DOUBLE n011_n111 = lerp(n011, n111, xs);

        DOUBLE n00x_n01x = lerp(n000_n100, n010_n110, ys);
        DOUBLE n10x_n11x = lerp(n001_n101, n011_n111, ys);

        return lerp(n00x_n01x, n10x_n11x, zs);
    }

    template<size_t VEC_LANES> inline typename fp2::simd::type_vec<double, VEC_LANES>::TYPE octaves3d(
            typename fp2::simd::type_vec<double, VEC_LANES>::TYPE x, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE y, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE z,
            typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE seed, size_t octaves) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;
        using DOUBLE_MASK = typename fp2::simd::type_vec<double, VEC_LANES>::BOOL;

        static const DOUBLE INITIAL_PERSISTENCE = 1.0;
        static const DOUBLE PERSISTENCE = 0.5;

        static const DOUBLE LANCULARITY = 2.0;

        DOUBLE value = 0.0;
        DOUBLE persistence = INITIAL_PERSISTENCE;

        for (size_t curOctave = 0; curOctave < octaves; curOctave++, persistence *= PERSISTENCE, x *= LANCULARITY, y *= LANCULARITY, z *= LANCULARITY) {
            DOUBLE nx = makeInt32Range<VEC_LANES>(x);
            DOUBLE ny = makeInt32Range<VEC_LANES>(y);
            DOUBLE nz = makeInt32Range<VEC_LANES>(z);

            value += noise3d<VEC_LANES>(nx, ny, nz, seed + curOctave) * persistence;
        }

        return value;
    }

    template<size_t VEC_LANES> inline typename fp2::simd::type_vec<double, VEC_LANES>::TYPE octaves3dVarying(
            typename fp2::simd::type_vec<double, VEC_LANES>::TYPE x, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE y, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE z,
            typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE seed, typename fp2::simd::type_vec<double, VEC_LANES>::UINT::TYPE octaves) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;
        using DOUBLE_MASK = typename fp2::simd::type_vec<double, VEC_LANES>::BOOL;
        using ULONG = typename fp2::simd::type_vec<double, VEC_LANES>::ULONG::TYPE;

        static const DOUBLE INITIAL_PERSISTENCE = 1.0;
        static const DOUBLE PERSISTENCE = 0.5;

        static const DOUBLE LANCULARITY = 2.0;

        DOUBLE value = 0.0;
        DOUBLE persistence = INITIAL_PERSISTENCE;

        for (size_t curOctave = 0, maxOctaves = horizontal_max(octaves); curOctave < maxOctaves; curOctave++, persistence *= PERSISTENCE, x *= LANCULARITY, y *= LANCULARITY, z *= LANCULARITY) {
            DOUBLE nx = makeInt32Range<VEC_LANES>(x);
            DOUBLE ny = makeInt32Range<VEC_LANES>(y);
            DOUBLE nz = makeInt32Range<VEC_LANES>(z);

            value = if_add(ULONG(curOctave) < extend(octaves), value, noise3d<VEC_LANES>(nx, ny, nz, seed + curOctave) * persistence);
        }

        return value;
    }

    template<size_t VEC_LANES> inline double octaves3dPoint(double _x, double _y, double _z, int32_t seed, size_t octaves) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;
        using DOUBLE_MASK = typename fp2::simd::type_vec<double, VEC_LANES>::BOOL;
        using INT = typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE;
        using ULONG = typename fp2::simd::type_vec<double, VEC_LANES>::ULONG::TYPE;

        static const DOUBLE INITIAL_LACUNARITY = fp2::simd::increment_shift<DOUBLE>();
        static const DOUBLE LANCULARITY = 1 << VEC_LANES;

        static const DOUBLE INITIAL_PERSISTENCE = 1.0 / fp2::simd::increment_shift<DOUBLE>();
        static const DOUBLE PERSISTENCE = 1.0 / (1 << VEC_LANES);

        static const INT INCREMENT = fp2::simd::increment<INT>();
        static const ULONG INCREMENT_L = fp2::simd::increment<ULONG>();

        DOUBLE x = _x * INITIAL_LACUNARITY;
        DOUBLE y = _y * INITIAL_LACUNARITY;
        DOUBLE z = _z * INITIAL_LACUNARITY;

        DOUBLE value = 0.0;
        DOUBLE persistence = INITIAL_PERSISTENCE;

        size_t curOctave = 0;
        for (; curOctave < (octaves & ~(VEC_LANES - 1)); curOctave += VEC_LANES, persistence *= PERSISTENCE, x *= LANCULARITY, y *= LANCULARITY, z *= LANCULARITY) {
            DOUBLE nx = makeInt32Range<VEC_LANES>(x);
            DOUBLE ny = makeInt32Range<VEC_LANES>(y);
            DOUBLE nz = makeInt32Range<VEC_LANES>(z);

            value += noise3d<VEC_LANES>(nx, ny, nz, seed + curOctave + INCREMENT) * persistence;
        }

        if (curOctave < octaves) { //there are some number of octaves left, let's do them with a mask
            DOUBLE nx = makeInt32Range<VEC_LANES>(x);
            DOUBLE ny = makeInt32Range<VEC_LANES>(y);
            DOUBLE nz = makeInt32Range<VEC_LANES>(z);

            value = if_add(curOctave + INCREMENT_L < octaves, value, noise3d<VEC_LANES>(nx, ny, nz, seed + curOctave + INCREMENT) * persistence);
        }

        return horizontal_add(value);
    }

    template<size_t VEC_LANES> inline typename fp2::simd::type_vec<double, VEC_LANES>::TYPE octaves2d(
            typename fp2::simd::type_vec<double, VEC_LANES>::TYPE x, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE z,
            typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE seed, size_t octaves) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;
        using DOUBLE_MASK = typename fp2::simd::type_vec<double, VEC_LANES>::BOOL;

        static const DOUBLE INITIAL_PERSISTENCE = 1.0;
        static const DOUBLE PERSISTENCE = 0.5;

        static const DOUBLE LANCULARITY = 2.0;

        DOUBLE value = 0.0;
        DOUBLE persistence = INITIAL_PERSISTENCE;

        for (size_t curOctave = 0; curOctave < octaves; curOctave++, persistence *= PERSISTENCE, x *= LANCULARITY, z *= LANCULARITY) {
            DOUBLE nx = makeInt32Range<VEC_LANES>(x);
            DOUBLE nz = makeInt32Range<VEC_LANES>(z);

            value += noise3d<VEC_LANES>(nx, 0.0, nz, seed + curOctave) * persistence;
        }

        return value;
    }

    template<size_t VEC_LANES> inline double octaves2dPoint(double _x, double _z, int32_t seed, size_t octaves) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;
        using DOUBLE_MASK = typename fp2::simd::type_vec<double, VEC_LANES>::BOOL;
        using INT = typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE;
        using ULONG = typename fp2::simd::type_vec<double, VEC_LANES>::ULONG::TYPE;

        static const DOUBLE INITIAL_LACUNARITY = fp2::simd::increment_shift<DOUBLE>();
        static const DOUBLE LANCULARITY = 1 << VEC_LANES;

        static const DOUBLE INITIAL_PERSISTENCE = 1.0 / fp2::simd::increment_shift<DOUBLE>();
        static const DOUBLE PERSISTENCE = 1.0 / (1 << VEC_LANES);

        static const INT INCREMENT = fp2::simd::increment<INT>();
        static const ULONG INCREMENT_L = fp2::simd::increment<ULONG>();

        DOUBLE x = _x * INITIAL_LACUNARITY;
        DOUBLE z = _z * INITIAL_LACUNARITY;

        DOUBLE value = 0.0;
        DOUBLE persistence = INITIAL_PERSISTENCE;

        size_t curOctave = 0;
        for (; curOctave < (octaves & ~(VEC_LANES - 1)); curOctave += VEC_LANES, persistence *= PERSISTENCE, x *= LANCULARITY, z *= LANCULARITY) {
            DOUBLE nx = makeInt32Range<VEC_LANES>(x);
            DOUBLE nz = makeInt32Range<VEC_LANES>(z);

            value += noise3d<VEC_LANES>(nx, 0.0, nz, seed + curOctave + INCREMENT) * persistence;
        }

        if (curOctave < octaves) { //there are some number of octaves left, let's do them with a mask
            DOUBLE nx = makeInt32Range<VEC_LANES>(x);
            DOUBLE nz = makeInt32Range<VEC_LANES>(z);

            value = if_add(curOctave + INCREMENT_L < octaves, value, noise3d<VEC_LANES>(nx, 0.0, nz, seed + curOctave + INCREMENT) * persistence);
        }

        return horizontal_add(value);
    }
}

#endif //CWG_NOISE_H
