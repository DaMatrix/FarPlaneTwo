/*
 * Adapted from The MIT License (MIT)
 *
 * Copyright (c) 2020-2020 DaPorkchop_
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

package net.daporkchop.fp2.strategy.heightmap;

import io.github.opencubicchunks.cubicchunks.cubicgen.customcubic.CustomCubicWorldType;
import io.netty.buffer.ByteBuf;
import lombok.AccessLevel;
import lombok.Getter;
import lombok.NonNull;
import net.daporkchop.fp2.FP2;
import net.daporkchop.fp2.FP2Config;
import net.daporkchop.fp2.strategy.common.IFarWorld;
import net.daporkchop.fp2.strategy.heightmap.gen.HeightmapGenerator;
import net.daporkchop.fp2.strategy.heightmap.gen.exact.CCHeightmapGenerator;
import net.daporkchop.fp2.strategy.heightmap.gen.exact.VanillaHeightmapGenerator;
import net.daporkchop.fp2.strategy.heightmap.gen.rough.CWGHeightmapGenerator;
import net.daporkchop.fp2.strategy.heightmap.scale.HeightmapScaler;
import net.daporkchop.fp2.strategy.heightmap.scale.HeightmapScalerMax;
import net.daporkchop.fp2.util.Constants;
import net.daporkchop.fp2.util.threading.CachedBlockAccess;
import net.daporkchop.lib.primitive.map.concurrent.ObjObjConcurrentHashMap;
import net.minecraft.world.WorldServer;

import java.util.Map;
import java.util.concurrent.CompletableFuture;

import static net.daporkchop.fp2.server.ServerConstants.*;
import static net.daporkchop.lib.common.util.PorkUtil.*;

/**
 * Representation of a world used by the heightmap rendering mode.
 *
 * @author DaPorkchop_
 */
@Getter
public class HeightmapWorld implements IFarWorld<HeightmapPos, HeightmapPiece> {
    protected final WorldServer world;

    protected final HeightmapGenerator generatorRough;
    protected final HeightmapGenerator generatorExact;
    protected final HeightmapScaler scaler;

    protected final HeightmapStorage storage;

    @Getter(AccessLevel.NONE)
    protected final Map<HeightmapPos, CompletableFuture<HeightmapPiece>> cache = new ObjObjConcurrentHashMap<>();

    public HeightmapWorld(@NonNull WorldServer world) {
        this.world = world;

        //TODO: this entire section should be re-done to utilize some form of registry so other mods can add their own rough generators

        //rough generator
        HeightmapGenerator roughGenerator = null;
        if (Constants.isCubicWorld(world)) {
            if (Constants.CWG && world.getWorldType() instanceof CustomCubicWorldType) {
                roughGenerator = new CWGHeightmapGenerator();
            }
        }

        //exact generator
        if (Constants.isCubicWorld(world)) {
            this.generatorExact = new CCHeightmapGenerator();
        } else {
            this.generatorExact = new VanillaHeightmapGenerator();
        }

        if (roughGenerator != null) { //rough generator has been set, so use it
            (this.generatorRough = roughGenerator).init(world);
        } else { //rough generator wasn't set, use the exact generator for this as well
            FP2.LOGGER.warn("No rough generator exists for world {} (type: {})! Falling back to exact generator, this will have serious performance implications.", world, world.getWorldType());
            this.generatorRough = this.generatorExact;
        }
        this.generatorExact.init(world);

        this.scaler = new HeightmapScalerMax();

        this.storage = new HeightmapStorage(world);
    }

    @Override
    public HeightmapPiece getPieceBlocking(@NonNull HeightmapPos pos) {
        return this.loadFullPiece(pos).join();
    }

    @Override
    public HeightmapPiece getPieceNowOrLoadAsync(@NonNull HeightmapPos pos) {
        return this.loadFullPiece(pos).getNow(null);
    }

    @Override
    public void blockChanged(int x, int y, int z) {
        //TODO: this whole thing is not really concurrency-safe
        /*this.loadFullPiece(key)
                .thenApplyAsync(piece -> {
                    if (this.dirtyChunks.replace(key, 1, 2)) {
                        CachedBlockAccess world = ((CachedBlockAccess.Holder) this.world).fp2_cachedBlockAccess();
                        try {
                            this.generator.generateExact(world, piece);
                            this.savePiece(piece);
                        } finally {
                            this.dirtyChunks.remove(key, 2);
                        }
                    }
                    return piece;
                }, GENERATION_WORKERS)
                .thenAcceptAsync(((IFarContext) this.world).fp2_tracker()::pieceChanged);*/
    }

    protected CompletableFuture<HeightmapPiece> loadFullPiece(HeightmapPos key) {
        return this.cache.computeIfAbsent(key, pos -> {
            CompletableFuture<HeightmapPiece> future = new CompletableFuture<>();
            GENERATION_WORKERS.submit(() -> {
                //load piece if possible
                HeightmapPiece piece = this.storage.loadAndUnpack(pos);
                if (piece != null) {
                    future.complete(piece);
                } else if (pos.level() == 0) { //root level, generate piece
                    try {
                        future.complete(this.generatePiece(pos.x(), pos.z(), this.generatorRough));
                    } catch (Exception e) {
                        future.completeExceptionally(e);
                    }
                } else {
                    CompletableFuture<HeightmapPiece>[] srcFutures = uncheckedCast(this.scaler.inputs(pos)
                            .map(this::loadFullPiece)
                            .toArray(CompletableFuture[]::new));

                    CompletableFuture.allOf(srcFutures)
                            .thenApplyAsync(v -> {
                                HeightmapPiece[] srcs = new HeightmapPiece[srcFutures.length];
                                for (int i = 0, len = srcs.length; i < len; i++) {
                                    srcs[i] = srcFutures[i].join();
                                }

                                HeightmapPiece dst = new HeightmapPiece(pos.x(), pos.z(), pos.level());

                                this.scaler.scale(srcs, dst);
                                this.savePiece(dst);
                                return dst;
                            }, GENERATION_WORKERS)
                            .whenComplete((v, t) -> {
                                if (t != null) {
                                    future.completeExceptionally(t);
                                } else {
                                    future.complete(v);
                                }
                            });
                }
            });
            return future;
        });
    }

    protected HeightmapPiece generatePiece(int x, int z, @NonNull HeightmapGenerator generator) {
        CachedBlockAccess world = ((CachedBlockAccess.Holder) this.world).cachedBlockAccess();
        HeightmapPiece piece = new HeightmapPiece(x, z, 0);
        generator.generate(world, piece);

        try {
            return piece;
        } finally {
            this.savePiece(piece);
        }
    }

    protected void savePiece(@NonNull HeightmapPiece piece) {
        if (FP2Config.debug.disableWrite || FP2Config.debug.disablePersistence || !piece.isDirty()) {
            return;
        }

        ByteBuf packed = null;

        try {
            piece.readLock().lock();
            try {
                packed = this.storage.pack(piece);
            } finally {
                piece.readLock().unlock();
            }

            this.storage.store(piece.pos(), packed);
        } finally {
            if (packed != null) {
                packed.release();
            }
        }
    }
}
