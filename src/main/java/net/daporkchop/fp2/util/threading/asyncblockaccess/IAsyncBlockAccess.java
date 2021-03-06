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

package net.daporkchop.fp2.util.threading.asyncblockaccess;

import lombok.NonNull;
import net.daporkchop.fp2.compat.vanilla.IBlockHeightAccess;
import net.daporkchop.fp2.util.threading.futurecache.GenerationNotAllowedException;
import net.minecraft.block.state.IBlockState;
import net.minecraft.tileentity.TileEntity;
import net.minecraft.util.EnumFacing;
import net.minecraft.util.math.BlockPos;
import net.minecraft.util.math.ChunkPos;
import net.minecraft.util.math.Vec3i;

import java.util.function.Function;
import java.util.stream.Stream;

/**
 * Thread-safe cache for block data in a world.
 *
 * @author DaPorkchop_
 */
public interface IAsyncBlockAccess extends IBlockHeightAccess {
    /**
     * Asynchronously prefetches the columns at the given positions into a single {@link IBlockHeightAccess}.
     * <p>
     * The returned {@link IBlockHeightAccess} will be able to respond to queries outside of the requested area, but they may be significantly
     * slower.
     *
     * @param columns a {@link Stream} containing the positions of all the columns to get
     * @return a single {@link IBlockHeightAccess} covering all of the given columns
     * @see #prefetchWithoutGenerating(Stream)
     */
    IBlockHeightAccess prefetch(@NonNull Stream<ChunkPos> columns);

    /**
     * Asynchronously prefetches the columns at the given positions into a single {@link IBlockHeightAccess}.
     * <p>
     * The returned {@link IBlockHeightAccess} will be able to respond to queries outside of the requested area, but they may be significantly
     * slower.
     *
     * @param columns a {@link Stream} containing the positions of all the columns to get
     * @return a single {@link IBlockHeightAccess} covering all of the given columns
     * @throws GenerationNotAllowedException if any terrain would need to be generated
     * @see #prefetch(Stream)
     */
    IBlockHeightAccess prefetchWithoutGenerating(@NonNull Stream<ChunkPos> columns) throws GenerationNotAllowedException;

    /**
     * Asynchronously prefetches the columns and cubes at the given positions into a single {@link IBlockHeightAccess}.
     * <p>
     * The returned {@link IBlockHeightAccess} will be able to respond to queries outside of the requested area, but they may be significantly
     * slower.
     *
     * @param columns              a {@link Stream} containing the positions of all the columns to get
     * @param cubesMappingFunction a function to produce a {@link Stream} containing the positions of all the cubes to get. The input parameter
     *                             is an {@link IBlockHeightAccess} containing only the requested chunks, but no cubes. Note that implementations
     *                             may choose to ignore this parameter.
     * @return a single {@link IBlockHeightAccess} covering all of the given columns and cubes
     * @see #prefetchWithoutGenerating(Stream, Function)
     */
    IBlockHeightAccess prefetch(@NonNull Stream<ChunkPos> columns, @NonNull Function<IBlockHeightAccess, Stream<Vec3i>> cubesMappingFunction);

    /**
     * Asynchronously prefetches the columns and cubes at the given positions into a single {@link IBlockHeightAccess}.
     * <p>
     * The returned {@link IBlockHeightAccess} will be able to respond to queries outside of the requested area, but they may be significantly
     * slower.
     *
     * @param columns              a {@link Stream} containing the positions of all the columns to get
     * @param cubesMappingFunction a function to produce a {@link Stream} containing the positions of all the cubes to get. The input parameter
     *                             is an {@link IBlockHeightAccess} containing only the requested chunks, but no cubes. Note that implementations
     *                             may choose to ignore this parameter.
     * @return a single {@link IBlockHeightAccess} covering all of the given columns and cubes
     * @throws GenerationNotAllowedException if any terrain would need to be generated
     * @see #prefetch(Stream, Function)
     */
    IBlockHeightAccess prefetchWithoutGenerating(@NonNull Stream<ChunkPos> columns, @NonNull Function<IBlockHeightAccess, Stream<Vec3i>> cubesMappingFunction) throws GenerationNotAllowedException;

    /**
     * @return whether or not any columns in the given tile exist
     */
    boolean anyColumnIntersects(int tileX, int tileZ, int level);

    /**
     * @return whether or not any cubes in the given area exist
     */
    boolean anyCubeIntersects(int tileX, int tileY, int tileZ, int level);

    @Override
    default boolean isAirBlock(BlockPos pos) {
        IBlockState state = this.getBlockState(pos);
        return state.getBlock().isAir(state, this, pos);
    }

    @Override
    default int getStrongPower(BlockPos pos, EnumFacing direction) {
        return this.getBlockState(pos).getStrongPower(this, pos, direction);
    }

    @Override
    default boolean isSideSolid(BlockPos pos, EnumFacing side, boolean _default) {
        return this.getBlockState(pos).isSideSolid(this, pos, side);
    }

    /**
     * @deprecated access to tile entities is not allowed
     */
    @Override
    @Deprecated
    default TileEntity getTileEntity(BlockPos pos) {
        throw new UnsupportedOperationException();
    }

    /**
     * Allows access to the {@link IAsyncBlockAccess} belonging to a {@link net.minecraft.world.WorldServer}.
     *
     * @author DaPorkchop_
     */
    interface Holder {
        IAsyncBlockAccess asyncBlockAccess();
    }
}
