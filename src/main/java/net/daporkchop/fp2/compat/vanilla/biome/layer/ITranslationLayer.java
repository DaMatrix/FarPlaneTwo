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

package net.daporkchop.fp2.compat.vanilla.biome.layer;

import lombok.NonNull;
import net.daporkchop.lib.common.pool.array.ArrayAllocator;

/**
 * A {@link IFastLayer} whose child requests are the same size as the initial input request.
 *
 * @author DaPorkchop_
 */
public interface ITranslationLayer extends IFastLayer {
    /**
     * @return the next layer in the generation chain
     */
    IFastLayer child();

    @Override
    default void getGrid(@NonNull ArrayAllocator<int[]> alloc, int x, int z, int sizeX, int sizeZ, @NonNull int[] out) {
        this.child().getGrid(alloc, x, z, sizeX, sizeZ, out);

        this.getGrid0(x, z, sizeX, sizeZ, out);
    }

    void getGrid0(int x, int z, int sizeX, int sizeZ, @NonNull int[] inout);

    @Override
    default void multiGetGrids(@NonNull ArrayAllocator<int[]> alloc, int x, int z, int size, int dist, int depth, int count, @NonNull int[] out) {
        this.child().multiGetGrids(alloc, x, z, size, dist, depth, count, out);

        this.multiGetGrids0(x, z, size, dist, depth, count, out);
    }

    void multiGetGrids0(int x, int z, int size, int dist, int depth, int count, @NonNull int[] inout);
}
