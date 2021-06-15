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

package net.daporkchop.fp2.compat.cc.asyncblockaccess;

import io.github.opencubicchunks.cubicchunks.api.util.CubePos;
import io.github.opencubicchunks.cubicchunks.api.world.IColumn;
import io.github.opencubicchunks.cubicchunks.api.world.ICube;
import io.github.opencubicchunks.cubicchunks.api.world.ICubeProviderServer;
import io.github.opencubicchunks.cubicchunks.api.world.ICubicWorldServer;
import io.github.opencubicchunks.cubicchunks.api.world.storage.ICubicStorage;
import io.github.opencubicchunks.cubicchunks.core.asm.mixin.ICubicWorldInternal;
import io.github.opencubicchunks.cubicchunks.core.server.chunkio.AsyncBatchingCubeIO;
import io.github.opencubicchunks.cubicchunks.core.server.chunkio.ICubeIO;
import io.github.opencubicchunks.cubicchunks.core.world.ICubeProviderInternal;
import io.github.opencubicchunks.cubicchunks.core.world.cube.Cube;
import lombok.NonNull;
import lombok.SneakyThrows;
import net.daporkchop.fp2.compat.cc.biome.Column2dBiomeAccessWrapper;
import net.daporkchop.fp2.compat.cc.biome.CubeBiomeAccessWrapper;
import net.daporkchop.fp2.compat.cc.cube.CubeWithoutWorld;
import net.daporkchop.fp2.compat.vanilla.IBiomeAccess;
import net.daporkchop.fp2.compat.vanilla.IBlockHeightAccess;
import net.daporkchop.fp2.server.worldlistener.IWorldChangeListener;
import net.daporkchop.fp2.server.worldlistener.WorldChangeListenerManager;
import net.daporkchop.fp2.util.datastructure.ConcurrentBooleanHashSegtreeInt;
import net.daporkchop.fp2.util.threading.ServerThreadExecutor;
import net.daporkchop.fp2.util.threading.asyncblockaccess.AsyncCacheNBTBase;
import net.daporkchop.fp2.util.threading.asyncblockaccess.IAsyncBlockAccess;
import net.daporkchop.fp2.util.threading.futurecache.GenerationNotAllowedException;
import net.daporkchop.fp2.util.threading.futurecache.IAsyncCache;
import net.daporkchop.fp2.util.threading.lazy.LazyFutureTask;
import net.daporkchop.lib.common.util.PorkUtil;
import net.daporkchop.lib.concurrent.PFutures;
import net.daporkchop.lib.unsafe.PUnsafe;
import net.minecraft.block.state.IBlockState;
import net.minecraft.nbt.NBTTagCompound;
import net.minecraft.util.math.BlockPos;
import net.minecraft.util.math.ChunkPos;
import net.minecraft.util.math.Vec3i;
import net.minecraft.world.EnumSkyBlock;
import net.minecraft.world.World;
import net.minecraft.world.WorldServer;
import net.minecraft.world.WorldType;
import net.minecraft.world.biome.Biome;
import net.minecraft.world.chunk.Chunk;
import net.minecraft.world.chunk.storage.ExtendedBlockStorage;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Stream;

import static net.daporkchop.lib.common.util.PorkUtil.*;

/**
 * Default implementation of {@link IAsyncBlockAccess} for vanilla worlds.
 *
 * @author DaPorkchop_
 */
public class CCAsyncBlockAccessImpl implements IAsyncBlockAccess, IWorldChangeListener {
    protected static final long ASYNCBATCHINGCUBEIO_STORAGE_OFFSET = PUnsafe.pork_getOffset(AsyncBatchingCubeIO.class, "storage");

    protected final WorldServer world;
    protected final ICubeIO io;
    protected final ICubicStorage storage;

    protected final ColumnCache columns = new ColumnCache();
    protected final CubeCache cubes = new CubeCache();

    protected final ExtendedBlockStorage emptyStorage;

    protected final ConcurrentBooleanHashSegtreeInt columnsExistCache;
    protected final ConcurrentBooleanHashSegtreeInt cubesExistCache;

    public CCAsyncBlockAccessImpl(@NonNull WorldServer world) {
        this.world = world;
        this.io = ((ICubeProviderInternal.Server) ((ICubicWorldInternal) world).getCubeCache()).getCubeIO();
        this.storage = PUnsafe.getObject(AsyncBatchingCubeIO.class.cast(this.io), ASYNCBATCHINGCUBEIO_STORAGE_OFFSET);

        this.emptyStorage = new ExtendedBlockStorage(0, world.provider.hasSkyLight());

        WorldChangeListenerManager.add(this.world, this);

        this.columnsExistCache = new ConcurrentBooleanHashSegtreeInt(2, () -> {
            List<int[]> positions = new ArrayList<>();
            CCAsyncBlockAccessImpl.this.storage.forEachColumn(pos -> positions.add(new int[]{ pos.x, pos.z }));
            return positions.stream();
        });
        this.cubesExistCache = new ConcurrentBooleanHashSegtreeInt(3, () -> {
            List<int[]> positions = new ArrayList<>();
            CCAsyncBlockAccessImpl.this.storage.forEachCube(pos -> positions.add(new int[]{ pos.getX(), pos.getY(), pos.getZ() }));
            return positions.stream();
        });
    }

    @Override
    public IBlockHeightAccess prefetch(@NonNull Stream<ChunkPos> columns) {
        //collect all futures into a list first in order to issue all tasks at once before blocking, thus ensuring maximum parallelism
        LazyFutureTask<IColumn>[] columnFutures = uncheckedCast(columns.map(pos -> this.columns.get(pos, true)).toArray(LazyFutureTask[]::new));

        return new PrefetchedColumnsCCAsyncBlockAccess(this, this.world, LazyFutureTask.scatterGather(columnFutures).stream());
    }

    @Override
    public IBlockHeightAccess prefetchWithoutGenerating(@NonNull Stream<ChunkPos> columns) throws GenerationNotAllowedException {
        //collect all futures into a list first in order to issue all tasks at once before blocking, thus ensuring maximum parallelism
        LazyFutureTask<IColumn>[] columnFutures = uncheckedCast(columns.map(pos -> this.columns.get(pos, false)).toArray(LazyFutureTask[]::new));

        return new PrefetchedColumnsCCAsyncBlockAccess(this, this.world, LazyFutureTask.scatterGather(columnFutures).stream()
                .peek(GenerationNotAllowedException.throwIfNull()));
    }

    @Override
    public IBlockHeightAccess prefetch(@NonNull Stream<ChunkPos> columns, @NonNull Function<IBlockHeightAccess, Stream<Vec3i>> cubesMappingFunction) {
        //collect all futures into a list first in order to issue all tasks at once before blocking, thus ensuring maximum parallelism
        LazyFutureTask<IColumn>[] columnFutures = uncheckedCast(columns.map(pos -> this.columns.get(pos, true)).toArray(LazyFutureTask[]::new));
        List<IColumn> columnList = LazyFutureTask.scatterGather(columnFutures);

        LazyFutureTask<ICube>[] cubeFutures = uncheckedCast(cubesMappingFunction.apply(new PrefetchedColumnsCCAsyncBlockAccess(this, this.world, columnList.stream()))
                .map(vec -> new CubePos(vec.getX(), vec.getY(), vec.getZ())).map(pos -> this.cubes.get(pos, true)).toArray(LazyFutureTask[]::new));

        return new PrefetchedCubesCCAsyncBlockAccess(this, this.world, columnList.stream(), LazyFutureTask.scatterGather(cubeFutures).stream());
    }

    @Override
    public IBlockHeightAccess prefetchWithoutGenerating(@NonNull Stream<ChunkPos> columns, @NonNull Function<IBlockHeightAccess, Stream<Vec3i>> cubesMappingFunction) throws GenerationNotAllowedException {
        //collect all futures into a list first in order to issue all tasks at once before blocking, thus ensuring maximum parallelism
        LazyFutureTask<IColumn>[] columnFutures = uncheckedCast(columns.map(pos -> this.columns.get(pos, false)).toArray(LazyFutureTask[]::new));
        List<IColumn> columnList = LazyFutureTask.scatterGather(columnFutures);
        columnList.forEach(GenerationNotAllowedException.throwIfNull());

        LazyFutureTask<ICube>[] cubeFutures = uncheckedCast(cubesMappingFunction.apply(new PrefetchedColumnsCCAsyncBlockAccess(this, this.world, columnList.stream()))
                .map(vec -> new CubePos(vec.getX(), vec.getY(), vec.getZ())).map(pos -> this.cubes.get(pos, false)).toArray(LazyFutureTask[]::new));

        return new PrefetchedCubesCCAsyncBlockAccess(this, this.world, columnList.stream(), LazyFutureTask.scatterGather(cubeFutures).stream()
                .peek(GenerationNotAllowedException.throwIfNull()));
    }

    @Override
    public void onColumnSaved(@NonNull World world, int columnX, int columnZ, @NonNull NBTTagCompound nbt) {
        this.columnsExistCache.set(columnX, columnZ);
        this.columns.notifyUpdate(new ChunkPos(columnX, columnZ), nbt);
    }

    @Override
    public void onCubeSaved(@NonNull World world, int cubeX, int cubeY, int cubeZ, @NonNull NBTTagCompound nbt) {
        this.cubesExistCache.set(cubeX, cubeY, cubeZ);
        this.cubes.notifyUpdate(new CubePos(cubeX, cubeY, cubeZ), nbt);
    }

    @Override
    public boolean anyColumnIntersects(int tileX, int tileZ, int level) {
        return this.columnsExistCache.isSet(level, tileX, tileZ);
    }

    @Override
    public boolean anyCubeIntersects(int tileX, int tileY, int tileZ, int level) {
        return this.cubesExistCache.isSet(level, tileX, tileY, tileZ);
    }

    protected IColumn getColumn(int columnX, int columnZ) {
        return this.columns.get(new ChunkPos(columnX, columnZ), true).join();
    }

    @Override
    public int getTopBlockY(int blockX, int blockZ) {
        return this.getColumn(blockX >> 4, blockZ >> 4).getOpacityIndex().getTopBlockY(blockX & 0xF, blockZ & 0xF);
    }

    @Override
    public int getTopBlockYBelow(int blockX, int blockY, int blockZ) {
        return this.getColumn(blockX >> 4, blockZ >> 4).getOpacityIndex().getTopBlockYBelow(blockX & 0xF, blockY, blockZ & 0xF);
    }

    protected ICube getCube(int cubeX, int cubeY, int cubeZ) {
        return this.cubes.get(new CubePos(cubeX, cubeY, cubeZ), true).join();
    }

    @Override
    public int getBlockLight(BlockPos pos) {
        if (!this.world.isValid(pos)) {
            return 0;
        } else {
            return this.getCube(pos.getX() >> 4, pos.getY() >> 4, pos.getZ() >> 4).getLightFor(EnumSkyBlock.BLOCK, pos);
        }
    }

    @Override
    public int getSkyLight(BlockPos pos) {
        if (!this.world.provider.hasSkyLight()) {
            return 0;
        } else if (!this.world.isValid(pos)) {
            return 15;
        } else {
            return this.getCube(pos.getX() >> 4, pos.getY() >> 4, pos.getZ() >> 4).getLightFor(EnumSkyBlock.SKY, pos);
        }
    }

    @Override
    public IBlockState getBlockState(BlockPos pos) {
        return this.getCube(pos.getX() >> 4, pos.getY() >> 4, pos.getZ() >> 4).getBlockState(pos);
    }

    @Override
    public Biome getBiome(BlockPos pos) {
        return this.getCube(pos.getX() >> 4, pos.getY() >> 4, pos.getZ() >> 4).getBiome(pos);
    }

    @Override
    public WorldType getWorldType() {
        return this.world.getWorldType();
    }

    /**
     * {@link IAsyncCache} for columns.
     *
     * @author DaPorkchop_
     */
    protected class ColumnCache extends AsyncCacheNBTBase<ChunkPos, Object, IColumn> {
        @Override
        protected IColumn parseNBT(@NonNull ChunkPos key, @NonNull Object param, @NonNull NBTTagCompound nbt) {
            ICubeIO.PartialData<Chunk> data = new ICubeIO.PartialData<>(null, nbt);
            CCAsyncBlockAccessImpl.this.io.loadColumnAsyncPart(data, key.x, key.z);
            return (IColumn) data.getObject();
        }

        @Override
        @SneakyThrows(IOException.class)
        protected IColumn loadFromDisk(@NonNull ChunkPos key, @NonNull Object param) {
            ICubeIO.PartialData<Chunk> data = CCAsyncBlockAccessImpl.this.io.loadColumnNbt(key.x, key.z);
            CCAsyncBlockAccessImpl.this.io.loadColumnAsyncPart(data, key.x, key.z);
            return (IColumn) data.getObject();
        }

        @Override
        protected void triggerGeneration(@NonNull ChunkPos key, @NonNull Object param) {
            //load and immediately save column on server thread
            PFutures.runAsync(() -> {
                Chunk column = ((ICubicWorldServer) CCAsyncBlockAccessImpl.this.world)
                        .getCubeCache().getColumn(key.x, key.z, ICubeProviderServer.Requirement.POPULATE);
                if (column != null && !column.isEmpty()) {
                    CCAsyncBlockAccessImpl.this.io.saveColumn(column);
                }
            }, ServerThreadExecutor.INSTANCE).join();
        }
    }

    /**
     * {@link IAsyncCache} for cubes.
     *
     * @author DaPorkchop_
     */
    protected class CubeCache extends AsyncCacheNBTBase<CubePos, Chunk, ICube> {
        @Override
        protected Chunk getParamFor(@NonNull CubePos key, boolean allowGeneration) {
            Object o = CCAsyncBlockAccessImpl.this.columns.get(key.chunkPos(), allowGeneration).join();
            return (Chunk) o;
        }

        @Override
        protected ICube parseNBT(@NonNull CubePos key, @NonNull Chunk param, @NonNull NBTTagCompound nbt) {
            ICubeIO.PartialData<ICube> data = new ICubeIO.PartialData<>(null, nbt);
            CCAsyncBlockAccessImpl.this.io.loadCubeAsyncPart(data, param, key.getY());
            return data.getObject() != null && data.getObject().isSurfaceTracked() ? data.getObject() : null;
        }

        @Override
        @SneakyThrows(IOException.class)
        protected ICube loadFromDisk(@NonNull CubePos key, @NonNull Chunk param) {
            ICubeIO.PartialData<ICube> data = CCAsyncBlockAccessImpl.this.io.loadCubeNbt(param, key.getY());
            CCAsyncBlockAccessImpl.this.io.loadCubeAsyncPart(data, param, key.getY());
            return data.getObject() != null && data.getObject().isSurfaceTracked() ? data.getObject() : null;
        }

        @Override
        protected void triggerGeneration(@NonNull CubePos key, @NonNull Chunk param) {
            PFutures.runAsync(() -> {
                //TODO: save column as well if needed
                ICube cube = ((ICubicWorldServer) CCAsyncBlockAccessImpl.this.world)
                        .getCubeCache().getCube(key.getX(), key.getY(), key.getZ(), ICubeProviderServer.Requirement.LIGHT);
                if (cube != null && cube.isSurfaceTracked()) {
                    CCAsyncBlockAccessImpl.this.io.saveCube((Cube) cube);
                }
            }, ServerThreadExecutor.INSTANCE).join();
        }

        @Override
        protected ICube bakeValue(@NonNull CubePos key, @NonNull Chunk param, @NonNull ICube value) {
            //override per-cube biomes if necessary
            IBiomeAccess biomeAccess = value instanceof Cube && ((Cube) value).getBiomeArray() != null
                    ? new CubeBiomeAccessWrapper(((Cube) value).getBiomeArray())
                    : new Column2dBiomeAccessWrapper(value.getColumn().getBiomeArray());

            return new CubeWithoutWorld(PorkUtil.fallbackIfNull(value.getStorage(), CCAsyncBlockAccessImpl.this.emptyStorage), biomeAccess, key);
        }
    }
}
