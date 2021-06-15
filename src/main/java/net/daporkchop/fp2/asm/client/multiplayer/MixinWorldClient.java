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

package net.daporkchop.fp2.asm.client.multiplayer;

import lombok.NonNull;
import net.daporkchop.fp2.mode.api.IFarPos;
import net.daporkchop.fp2.mode.api.IFarRenderMode;
import net.daporkchop.fp2.mode.api.ctx.IFarClientContext;
import net.daporkchop.fp2.mode.api.ctx.IFarWorldClient;
import net.daporkchop.fp2.mode.api.IFarTile;
import net.daporkchop.fp2.util.threading.ClientThreadExecutor;
import net.daporkchop.lib.primitive.map.concurrent.ObjObjConcurrentHashMap;
import net.minecraft.client.multiplayer.WorldClient;
import net.minecraft.world.World;
import org.spongepowered.asm.mixin.Implements;
import org.spongepowered.asm.mixin.Interface;
import org.spongepowered.asm.mixin.Mixin;
import org.spongepowered.asm.mixin.Unique;

import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.function.Function;

import static net.daporkchop.fp2.util.Constants.*;
import static net.daporkchop.lib.common.util.PorkUtil.*;

/**
 * @author DaPorkchop_
 */
@Mixin(WorldClient.class)
@Implements({
        @Interface(iface = IFarWorldClient.class, prefix = "fp2_world$", unique = true)
})
public abstract class MixinWorldClient extends World implements IFarWorldClient {
    @Unique
    private final Map<IFarRenderMode, IFarClientContext> contexts = new ObjObjConcurrentHashMap<>();
    @Unique
    private final Function<IFarRenderMode, IFarClientContext> computeFunction = m ->
            //always create context on client thread
            CompletableFuture.supplyAsync(() -> m.clientContext(uncheckedCast(this)), ClientThreadExecutor.INSTANCE).join();

    @Unique
    private IFarClientContext active;

    protected MixinWorldClient() {
        super(null, null, null, null, false);
    }

    @Override
    public <POS extends IFarPos, T extends IFarTile> IFarClientContext<POS, T> contextFor(@NonNull IFarRenderMode<POS, T> mode) {
        return uncheckedCast(this.contexts.computeIfAbsent(mode, this.computeFunction));
    }

    @Override
    public void switchTo(IFarRenderMode<?, ?> mode) {
        FP2_LOG.info("switching render mode to {}", mode != null ? mode.name() : null);
        this.active = mode != null ? this.contextFor(mode) : null;
    }

    @Override
    public <POS extends IFarPos, T extends IFarTile> IFarClientContext<POS, T> activeContext() {
        return uncheckedCast(this.active);
    }

    @Override
    public void close() {
        this.contexts.forEach((mode, context) -> context.close());
        this.contexts.clear();
    }
}
