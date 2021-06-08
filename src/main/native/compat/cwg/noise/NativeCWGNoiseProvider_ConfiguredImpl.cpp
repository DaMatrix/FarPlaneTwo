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

#include "cwg_noise.h"

class state_t {
public:
    double _heightVariationFactor;
    double _specialHeightVariationFactorBelowAverageY;
    double _heightVariationOffset;
    double _heightFactor;
    double _heightOffset;

    int32_t _selectorNoiseSeed;
    uint32_t _selectorNoiseOctaves;
    double _selectorNoiseFrequencyX;
    double _selectorNoiseFrequencyY;
    double _selectorNoiseFrequencyZ;
    double _selectorNoiseFactor;
    double _selectorNoiseOffset;

    int32_t _lowNoiseSeed;
    uint32_t _lowNoiseOctaves;
    double _lowNoiseFrequencyX;
    double _lowNoiseFrequencyY;
    double _lowNoiseFrequencyZ;
    double _lowNoiseFactor;
    double _lowNoiseOffset;

    int32_t _highNoiseSeed;
    uint32_t _highNoiseOctaves;
    double _highNoiseFrequencyX;
    double _highNoiseFrequencyY;
    double _highNoiseFrequencyZ;
    double _highNoiseFactor;
    double _highNoiseOffset;

    int32_t _depthNoiseSeed;
    uint32_t _depthNoiseOctaves;
    double _depthNoiseFrequencyX;
    double _depthNoiseFrequencyZ;
    double _depthNoiseFactor;
    double _depthNoiseOffset;

    fp2::simd::type_vec<int32_t, 4>::TYPE _allNoiseSeeds;
    fp2::simd::type_vec<uint32_t, 4>::TYPE _allNoiseOctaves;
    fp2::simd::type_vec<uint32_t, 4>::TYPE _allNoiseOctavesNoDepth;
    fp2::simd::type_vec<double, 4>::TYPE _allNoiseFrequenciesX;
    fp2::simd::type_vec<double, 4>::TYPE _allNoiseFrequenciesY;
    fp2::simd::type_vec<double, 4>::TYPE _allNoiseFrequenciesZ;
    fp2::simd::type_vec<double, 4>::TYPE _allNoiseFactors;
    fp2::simd::type_vec<double, 4>::TYPE _allNoiseOffsets;

    inline void setupVectorFields() {
        _allNoiseSeeds = fp2::simd::type_vec<int32_t, 4>::TYPE(_selectorNoiseSeed, _lowNoiseSeed, _highNoiseSeed, _depthNoiseSeed);
        _allNoiseOctaves = fp2::simd::type_vec<uint32_t, 4>::TYPE(_selectorNoiseOctaves, _lowNoiseOctaves, _highNoiseOctaves, _depthNoiseOctaves);
        _allNoiseOctavesNoDepth = fp2::simd::type_vec<uint32_t, 4>::TYPE(_selectorNoiseOctaves, _lowNoiseOctaves, _highNoiseOctaves, 0);
        _allNoiseFrequenciesX = fp2::simd::type_vec<double, 4>::TYPE(_selectorNoiseFrequencyX, _lowNoiseFrequencyX, _highNoiseFrequencyX, _depthNoiseFrequencyX);
        _allNoiseFrequenciesY = fp2::simd::type_vec<double, 4>::TYPE(_selectorNoiseFrequencyY, _lowNoiseFrequencyY, _highNoiseFrequencyY, 0.0);
        _allNoiseFrequenciesZ = fp2::simd::type_vec<double, 4>::TYPE(_selectorNoiseFrequencyZ, _lowNoiseFrequencyZ, _highNoiseFrequencyZ, _depthNoiseFrequencyZ);
        _allNoiseFactors = fp2::simd::type_vec<double, 4>::TYPE(_selectorNoiseFactor, _lowNoiseFactor, _highNoiseFactor, _depthNoiseFactor);
        _allNoiseOffsets = fp2::simd::type_vec<double, 4>::TYPE(_selectorNoiseOffset, _lowNoiseOffset, _highNoiseOffset, _depthNoiseOffset);
    }

    template<size_t VEC_LANES, bool USE_DEPTH_ARG> inline double generateAndMixAllNoise(double height, double variation, double depth, int32_t x, int32_t y, int32_t z) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;

        union {
            struct {
                double selector;
                double low;
                double high;
                double depth;
            } scalar;
            double vector[4];
        } noise;

        if constexpr (VEC_LANES == 4) {
            //there are 3 or 4 values to be computed, and we have 4 vector lanes available: perfect!
            if constexpr (USE_DEPTH_ARG) {
                //compute all 4 noise types at once, with depth noise being some irrelevant value
                DOUBLE noise_vec = fp2::cwg::noise::octaves3dVarying<VEC_LANES>(x * _allNoiseFrequenciesX, y * _allNoiseFrequenciesY, z * _allNoiseFrequenciesZ, _allNoiseSeeds, _allNoiseOctavesNoDepth);
                noise_vec = noise_vec * _allNoiseFactors + _allNoiseOffsets;
                noise_vec.store(&noise.vector[0]);

                //override depth value with input value
                noise.scalar.depth = depth;
            } else {
                //compute all 4 noise types at once
                DOUBLE noise_vec = fp2::cwg::noise::octaves3dVarying<VEC_LANES>(x * _allNoiseFrequenciesX, y * _allNoiseFrequenciesY, z * _allNoiseFrequenciesZ, _allNoiseSeeds, _allNoiseOctaves);
                noise_vec = noise_vec * _allNoiseFactors + _allNoiseOffsets;
                noise_vec.store(&noise.vector[0]);

                //post-process depth
                noise.scalar.depth = processDepthNoise<false>(noise.scalar.depth);
            }
        } else {
            //there are 3 or 4 values to be computed, but there are more than 4 lanes available: don't waste lanes by computing more values than we need to, let's vectorize on the octave count instead
            noise.scalar.selector = fp2::cwg::noise::octaves3dPoint<VEC_LANES>(x * _selectorNoiseFrequencyX, y * _selectorNoiseFrequencyY, z * _selectorNoiseFrequencyZ, _selectorNoiseSeed, _selectorNoiseOctaves)
                    * _selectorNoiseFactor + _selectorNoiseOffset;
            noise.scalar.low = fp2::cwg::noise::octaves3dPoint<VEC_LANES>(x * _lowNoiseFrequencyX, y * _lowNoiseFrequencyY, z * _lowNoiseFrequencyZ, _lowNoiseSeed, _lowNoiseOctaves)
                    * _lowNoiseFactor + _lowNoiseOffset;
            noise.scalar.high = fp2::cwg::noise::octaves3dPoint<VEC_LANES>(x * _highNoiseFrequencyX, y * _highNoiseFrequencyY, z * _highNoiseFrequencyZ, _highNoiseSeed, _highNoiseOctaves)
                    * _highNoiseFactor + _highNoiseOffset;

            if constexpr (USE_DEPTH_ARG) {
                //use input depth
                noise.scalar.depth = depth;
            } else {
                //compute depth from input
                noise.scalar.depth = fp2::cwg::noise::octaves2dPoint<VEC_LANES>(x * _depthNoiseFrequencyX, z * _depthNoiseFrequencyZ, _depthNoiseSeed, _depthNoiseOctaves)
                        * _depthNoiseFactor + _depthNoiseOffset;

                //post-process depth
                noise.scalar.depth = processDepthNoise<false>(noise.scalar.depth);
            }
        }

        height = height * _heightFactor + _heightOffset;
        variation *= _heightVariationFactor;
        if (height > y) {
            variation *= _specialHeightVariationFactorBelowAverageY;
        }
        variation += _heightVariationOffset;

        double d = fp2::cwg::noise::lerp(noise.scalar.low, noise.scalar.high, fp2::cwg::noise::clamp(noise.scalar.selector, 0.0, 1.0)) + noise.scalar.depth;
        d = d * variation + height;
        return d - fp2::cwg::noise::signum(variation) * y;
    }

    template<size_t VEC_LANES, bool USE_DEPTH_ARG> inline typename fp2::simd::type_vec<double, VEC_LANES>::TYPE generateAndMixAllNoise(
            typename fp2::simd::type_vec<double, VEC_LANES>::TYPE height, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE variation, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE depth,
            typename fp2::simd::type_vec<double, VEC_LANES>::TYPE x, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE y, typename fp2::simd::type_vec<double, VEC_LANES>::TYPE z) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;

        DOUBLE selector = fp2::cwg::noise::octaves3d<VEC_LANES>(x * _selectorNoiseFrequencyX, y * _selectorNoiseFrequencyY, z * _selectorNoiseFrequencyZ, _selectorNoiseSeed, _selectorNoiseOctaves)
                * _selectorNoiseFactor + _selectorNoiseOffset;
        DOUBLE low = fp2::cwg::noise::octaves3d<VEC_LANES>(x * _lowNoiseFrequencyX, y * _lowNoiseFrequencyY, z * _lowNoiseFrequencyZ, _lowNoiseSeed, _lowNoiseOctaves)
                * _lowNoiseFactor + _lowNoiseOffset;
        DOUBLE high = fp2::cwg::noise::octaves3d<VEC_LANES>(x * _highNoiseFrequencyX, y * _highNoiseFrequencyY, z * _highNoiseFrequencyZ, _highNoiseSeed, _highNoiseOctaves)
                * _highNoiseFactor + _highNoiseOffset;

        if constexpr (!USE_DEPTH_ARG) {
            depth = processDepthNoise<VEC_LANES>(fp2::cwg::noise::octaves2d<VEC_LANES>(x * _depthNoiseFrequencyX, z * _depthNoiseFrequencyZ, _depthNoiseSeed, _depthNoiseOctaves));
        }

        height = height * _heightFactor + _heightOffset;
        variation = if_mul(height > y, variation * _heightVariationFactor, _specialHeightVariationFactorBelowAverageY) + _heightVariationOffset;

        DOUBLE d = fp2::cwg::noise::lerp(low, high, max(min(selector, 1.0), 0.0)) + depth;
        d = d * variation + height;
        return d - select(variation == 0.0, 0.0, sign_combine(1.0, variation)) * y;
    }

    template<bool SCALE = true> inline double processDepthNoise(double depth) {
        if constexpr (SCALE) {
            depth = depth * _depthNoiseFactor + _depthNoiseOffset;
        }
        depth *= depth < 0.0 ? -0.9 : 3.0;
        depth -= 2.0;
        depth = fp2::cwg::noise::clamp(depth * (depth < 0.0 ? 5.0 / 28.0 : 0.125), -5.0 / 14.0, 0.125);
        depth *= 0.2 * 17.0 / 64.0;
        return depth;
    }

    template<size_t VEC_LANES> inline typename fp2::simd::type_vec<double, VEC_LANES>::TYPE processDepthNoise(typename fp2::simd::type_vec<double, VEC_LANES>::TYPE depth) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;

        depth = depth * _depthNoiseFactor + _depthNoiseOffset;
        depth *= select(depth < 0.0, DOUBLE(-0.9), DOUBLE(3.0));
        depth -= 2.0;
        depth = min(max(depth * select(depth < 0.0, DOUBLE(5.0 / 28.0), DOUBLE(0.125)), -5.0 / 14.0), 0.125);
        depth *= 0.2 * 17.0 / 64.0;
        return depth;
    }

    inline double generateDepthSingle(int32_t x, int32_t z) {
        constexpr size_t VEC_LANES = fp2::simd::LANES_32AND64;

        return processDepthNoise(fp2::cwg::noise::octaves2dPoint<VEC_LANES>(x * _depthNoiseFrequencyX, z * _depthNoiseFrequencyZ, _depthNoiseSeed, _depthNoiseOctaves));
    }

    inline void generateDepth2d(fp2::pinned_double_array& out, int32_t baseX, int32_t baseZ, int32_t scaleX, int32_t scaleZ, int32_t sizeX, int32_t sizeZ) {
        constexpr size_t VEC_LANES = fp2::simd::LANES_32AND64;

        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;
        using DOUBLE_MASK = typename fp2::simd::type_vec<double, VEC_LANES>::BOOL;
        using INT = typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE;
        using INT_MASK = typename fp2::simd::type_vec<double, VEC_LANES>::INT::BOOL;

        if ((sizeZ < VEC_LANES && sizeX != 1)) { //if sizeZ is less than the number of vector lanes, we can't do vectorized multi-dimensional iteration in SIMD
            for (int32_t i = 0, dx = 0; dx < sizeX; dx++) {
                for (int32_t dz = 0; dz < sizeZ; dz++, i++) {
                    out[i] = processDepthNoise(fp2::cwg::noise::octaves2dPoint<VEC_LANES>((baseX + dx * scaleX) * _depthNoiseFrequencyX, (baseZ + dz * scaleZ) * _depthNoiseFrequencyX, _depthNoiseSeed, _depthNoiseOctaves));
                }
            }
        } else {
            static const INT INCREMENT = fp2::simd::increment<INT>();

            INT x = baseX;
            INT z = baseZ + INCREMENT * scaleZ;

            const INT maxZ = baseZ + sizeZ * scaleZ;
            const INT stepZ = (int32_t) VEC_LANES * scaleZ;
            const INT resetZ = sizeZ * scaleZ;
            const INT stepX = scaleX;

            const size_t totalCount = sizeX * sizeZ;
            size_t index = 0;

            for (; index < (totalCount & ~(VEC_LANES - 1)); index += VEC_LANES) {
                processDepthNoise<VEC_LANES>(fp2::cwg::noise::octaves2d<VEC_LANES>(to_double(x) * _depthNoiseFrequencyX, to_double(z) * _depthNoiseFrequencyZ, _depthNoiseSeed, _depthNoiseOctaves))
                        .store(&out[index]);

                //increment z coordinates, resetting them and incrementing x if they reach the maximum value
                z += stepZ;
                INT_MASK ge = z >= maxZ;
                z = if_sub(ge, z, resetZ);
                x = if_add(ge, x, stepX);
            }

            if (index < totalCount) { //the number of samples remaining are less than the number of vector lanes, let's finish 'em up
                //x and z are already set up correctly
                processDepthNoise<VEC_LANES>(fp2::cwg::noise::octaves2d<VEC_LANES>(to_double(x) * _depthNoiseFrequencyX, to_double(z) * _depthNoiseFrequencyZ, _depthNoiseSeed, _depthNoiseOctaves))
                        .store_partial(totalCount & (VEC_LANES - 1), &out[index]);
            }
        }
    }

    template<size_t VEC_LANES> inline typename fp2::simd::type_vec<double, VEC_LANES>::TYPE loadWithWrap(
            fp2::pinned_double_array& arr, typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE i) {
        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;
        using INT = typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE;
        using INT_BOOL = typename fp2::simd::type_vec<double, VEC_LANES>::INT::BOOL;

        //using intrinsics here because vectorclass doesn't have a way to lookup 64-bit values with 32-bit indices
        if constexpr (INSTRSET >= 8 && VEC_LANES == 4) { //AVX2
            return DOUBLE(_mm256_i32gather_pd(&arr[0], i, 8));
        } else if constexpr (INSTRSET >= 9 && VEC_LANES == 8) { //AVX512
            return DOUBLE(_mm512_i32gather_pd(i, &arr[0], 8));
        } else { //AVX or SSE
            static const INT INCREMENT = fp2::simd::increment<INT>();

            int32_t _i[VEC_LANES];
            i.store(&_i[0]);

            INT_BOOL equals = (i - _i[0]) == INCREMENT;
            if (horizontal_and(equals)) { //indices are sequential
                return DOUBLE().load(&arr[_i[0]]);
            } else { //slow non-sequential lookup
                return lookup<(1 << 30)>(extend(i), &arr[0]);
            }
        }
    }

    template<bool USE_DEPTH_ARG> inline void generate3d(
            fp2::pinned_double_array& height, fp2::pinned_double_array& variation, fp2::pinned_double_array& depth, fp2::pinned_double_array& out,
            int32_t baseX, int32_t baseY, int32_t baseZ, int32_t scaleX, int32_t scaleY, int32_t scaleZ, int32_t sizeX, int32_t sizeY, int32_t sizeZ) {
        constexpr size_t VEC_LANES = fp2::simd::LANES_32AND64;

        using DOUBLE = typename fp2::simd::type_vec<double, VEC_LANES>::TYPE;
        using DOUBLE_MASK = typename fp2::simd::type_vec<double, VEC_LANES>::BOOL;
        using INT = typename fp2::simd::type_vec<double, VEC_LANES>::INT::TYPE;
        using INT_MASK = typename fp2::simd::type_vec<double, VEC_LANES>::INT::BOOL;

        if (sizeZ < VEC_LANES && (sizeY != 1 || sizeX != 1)) { //if sizeZ is less than the number of vector lanes, we can't do vectorized multi-dimensional iteration in SIMD
            for (int32_t i = 0, dx = 0; dx < sizeX; dx++) {
                for (int32_t dy = 0; dy < sizeY; dy++) {
                    for (int32_t inIdx = dx * sizeZ, dz = 0; dz < sizeZ; dz++, inIdx++, i++) {
                        out[i] = generateAndMixAllNoise<VEC_LANES, USE_DEPTH_ARG>(
                                height[inIdx], variation[inIdx], USE_DEPTH_ARG ? depth[inIdx] : 0.0,
                                baseX + dx * scaleX, baseY + dy * scaleY, baseZ + dz * scaleZ);
                    }
                }
            }
        } else {
            static const INT INCREMENT = fp2::simd::increment<INT>();

            INT x = baseX;
            INT y = baseY;
            INT z = baseZ + INCREMENT * scaleZ;
            INT readIndex = INCREMENT;

            const INT maxZ = baseZ + sizeZ * scaleZ;
            const INT stepZ = (int32_t) VEC_LANES * scaleZ;
            const INT resetZ = sizeZ * scaleZ;
            const INT stepY = scaleY;
            const INT maxY = baseY + sizeY * scaleY;
            const INT resetY = sizeY * scaleY;
            const INT stepX = scaleX;
            const INT stepReadIndex = (int32_t) VEC_LANES;

            const size_t totalCount = sizeX * sizeY * sizeZ;
            size_t index = 0;

            for (; index < (totalCount & ~(VEC_LANES - 1)); index += VEC_LANES) {
                generateAndMixAllNoise<VEC_LANES, USE_DEPTH_ARG>(
                        loadWithWrap<VEC_LANES>(height, readIndex),
                        loadWithWrap<VEC_LANES>(variation, readIndex),
                        USE_DEPTH_ARG ? loadWithWrap<VEC_LANES>(depth, readIndex) : 0.0,
                        to_double(x), to_double(y), to_double(z)).store(&out[index]);

                //increment z coordinates, resetting them and incrementing y if they reach the maximum value
                z += stepZ;
                readIndex += stepReadIndex;
                INT_MASK ge = z >= maxZ;
                z = if_sub(ge, z, resetZ);
                readIndex = if_sub(ge, readIndex, sizeZ);

                //increment y coordinates, resetting them and incrementing x if they reach the maximum value
                y = if_add(ge, y, stepY);
                ge = y >= maxY;
                y = if_sub(ge, y, resetY);
                x = if_add(ge, x, stepX);
                readIndex = if_add(ge, readIndex, sizeZ);
            }

            if (index < totalCount) { //the number of samples remaining are less than the number of vector lanes, let's finish 'em up
                //x and z are already set up correctly
                size_t remaining = totalCount & (VEC_LANES - 1);
                size_t baseIndex = readIndex[0];
                generateAndMixAllNoise<VEC_LANES, USE_DEPTH_ARG>(
                        DOUBLE().load_partial(remaining, &height[baseIndex]),
                        DOUBLE().load_partial(remaining, &variation[baseIndex]),
                        USE_DEPTH_ARG ? DOUBLE().load_partial(remaining, &depth[baseIndex]) : 0.0,
                        to_double(x), to_double(y), to_double(z)).store_partial(remaining, &out[index]);
            }
        }
    }

    inline double generateSingle(double height, double variation, int32_t x, int32_t y, int32_t z) {
        return generateAndMixAllNoise<fp2::simd::LANES_32AND64, false>(height, variation, 0.0, x, y, z);
    }

    inline double generateSingle(double height, double variation, double depth, int32_t x, int32_t y, int32_t z) {
        return generateAndMixAllNoise<fp2::simd::LANES_32AND64, true>(height, variation, depth, x, y, z);
    }
};

FP2_JNI(jlong, NativeCWGNoiseProvider_00024ConfiguredImpl, createState0)(JNIEnv* env, jclass cla,
        jdouble heightVariationFactor, jdouble specialHeightVariationFactorBelowAverageY, jdouble heightVariationOffset, jdouble heightFactor, jdouble heightOffset,
        jdouble selectorNoiseFactor, jdouble selectorNoiseOffset, jdouble selectorNoiseFrequencyX, jdouble selectorNoiseFrequencyY, jdouble selectorNoiseFrequencyZ, jint selectorNoiseSeed, jint selectorNoiseOctaves, jdouble selectorNoiseScale,
        jdouble lowNoiseFactor, jdouble lowNoiseOffset, jdouble lowNoiseFrequencyX, jdouble lowNoiseFrequencyY, jdouble lowNoiseFrequencyZ, jint lowNoiseSeed, jint lowNoiseOctaves, jdouble lowNoiseScale,
        jdouble highNoiseFactor, jdouble highNoiseOffset, jdouble highNoiseFrequencyX, jdouble highNoiseFrequencyY, jdouble highNoiseFrequencyZ, jint highNoiseSeed, jint highNoiseOctaves, jdouble highNoiseScale,
        jdouble depthNoiseFactor, jdouble depthNoiseOffset, jdouble depthNoiseFrequencyX, jdouble depthNoiseFrequencyZ, jint depthNoiseSeed, jint depthNoiseOctaves, jdouble depthNoiseScale) {
    state_t* state = new state_t;

    state->_heightVariationFactor = heightVariationFactor;
    state->_specialHeightVariationFactorBelowAverageY = specialHeightVariationFactorBelowAverageY;
    state->_heightVariationOffset = heightVariationOffset;
    state->_heightFactor = heightFactor;
    state->_heightOffset = heightOffset;

    state->_depthNoiseFactor = depthNoiseScale * depthNoiseFactor;
    state->_depthNoiseOffset = depthNoiseOffset - depthNoiseFactor;
    state->_depthNoiseFrequencyX = depthNoiseFrequencyX;
    state->_depthNoiseFrequencyZ = depthNoiseFrequencyZ;
    state->_depthNoiseSeed = depthNoiseSeed;
    state->_depthNoiseOctaves = (uint32_t) depthNoiseOctaves;

    state->_selectorNoiseFactor = selectorNoiseScale * selectorNoiseFactor;
    state->_selectorNoiseOffset = selectorNoiseOffset - selectorNoiseFactor;
    state->_selectorNoiseFrequencyX = selectorNoiseFrequencyX;
    state->_selectorNoiseFrequencyY = selectorNoiseFrequencyY;
    state->_selectorNoiseFrequencyZ = selectorNoiseFrequencyZ;
    state->_selectorNoiseSeed = selectorNoiseSeed;
    state->_selectorNoiseOctaves = (uint32_t) selectorNoiseOctaves;

    state->_lowNoiseFactor = lowNoiseScale * lowNoiseFactor;
    state->_lowNoiseOffset = lowNoiseOffset - lowNoiseFactor;
    state->_lowNoiseFrequencyX = lowNoiseFrequencyX;
    state->_lowNoiseFrequencyY = lowNoiseFrequencyY;
    state->_lowNoiseFrequencyZ = lowNoiseFrequencyZ;
    state->_lowNoiseSeed = lowNoiseSeed;
    state->_lowNoiseOctaves = (uint32_t) lowNoiseOctaves;

    state->_highNoiseFactor = highNoiseScale * highNoiseFactor;
    state->_highNoiseOffset = highNoiseOffset - highNoiseFactor;
    state->_highNoiseFrequencyX = highNoiseFrequencyX;
    state->_highNoiseFrequencyY = highNoiseFrequencyY;
    state->_highNoiseFrequencyZ = highNoiseFrequencyZ;
    state->_highNoiseSeed = highNoiseSeed;
    state->_highNoiseOctaves = (uint32_t) highNoiseOctaves;

    state->setupVectorFields();

    return (jlong) state;
}

FP2_JNI(void, NativeCWGNoiseProvider_00024ConfiguredImpl, deleteState0)(JNIEnv* env, jclass cla,
        jlong _state) {
    delete (state_t*) _state;
}

FP2_JNI(void, NativeCWGNoiseProvider_00024ConfiguredImpl, generateDepth2d0)(JNIEnv* env, jobject obj,
        jdoubleArray _out, jint baseX, jint baseZ, jint scaleX, jint scaleZ, jint sizeX, jint sizeZ, jlong _state) {
    fp2::pinned_double_array out(env, _out);
    ((state_t*) _state)->generateDepth2d(out, baseX, baseZ, scaleX, scaleZ, sizeX, sizeZ);
}

FP2_JNI(jdouble, NativeCWGNoiseProvider_00024ConfiguredImpl, generateDepthSingle0)(JNIEnv* env, jobject obj,
        jint x, jint z, jlong _state) {
    return ((state_t*) _state)->generateDepthSingle(x, z);
}

FP2_JNI(jdouble, NativeCWGNoiseProvider_00024ConfiguredImpl, generateSingle0noDepth)(JNIEnv* env, jobject obj,
        jdouble height, jdouble variation, jint x, jint y, jint z, jlong _state) {
    return ((state_t*) _state)->generateSingle(height, variation, x, y, z);
}

FP2_JNI(void, NativeCWGNoiseProvider_00024ConfiguredImpl, generate3d0noDepth)(JNIEnv* env, jobject obj,
        jdoubleArray _height, jdoubleArray _variation, jdoubleArray _out, jint baseX, jint baseY, jint baseZ, jint scaleX, jint scaleY, jint scaleZ, jint sizeX, jint sizeY, jint sizeZ, jlong _state) {
    fp2::pinned_double_array height(env, _height);
    fp2::pinned_double_array variation(env, _variation);
    fp2::pinned_double_array out(env, _out);
    ((state_t*) _state)->generate3d<false>(height, variation, variation, out, baseX, baseY, baseZ, scaleX, scaleY, scaleZ, sizeX, sizeY, sizeZ);
}

FP2_JNI(void, NativeCWGNoiseProvider_00024ConfiguredImpl, generate3d0depth)(JNIEnv* env, jobject obj,
        jdoubleArray _height, jdoubleArray _variation, jdoubleArray _depth, jdoubleArray _out, jint baseX, jint baseY, jint baseZ, jint scaleX, jint scaleY, jint scaleZ, jint sizeX, jint sizeY, jint sizeZ, jlong _state) {
    fp2::pinned_double_array height(env, _height);
    fp2::pinned_double_array variation(env, _variation);
    fp2::pinned_double_array depth(env, _depth);
    fp2::pinned_double_array out(env, _out);
    ((state_t*) _state)->generate3d<true>(height, variation, depth, out, baseX, baseY, baseZ, scaleX, scaleY, scaleZ, sizeX, sizeY, sizeZ);
}

FP2_JNI(jdouble, NativeCWGNoiseProvider_00024ConfiguredImpl, generateSingle0depth)(JNIEnv* env, jobject obj,
        jdouble height, jdouble variation, jdouble depth, jint x, jint y, jint z, jlong _state) {
    return ((state_t*) _state)->generateSingle(height, variation, depth, x, y, z);
}
