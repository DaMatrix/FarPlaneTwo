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

package net.daporkchop.fp2.mode.common.client.strategy;

import lombok.Getter;
import net.daporkchop.fp2.client.gl.commandbuffer.IDrawCommandBuffer;
import net.daporkchop.fp2.mode.api.IFarPos;
import net.daporkchop.fp2.mode.api.IFarTile;

import java.util.stream.IntStream;

import static net.daporkchop.fp2.mode.common.client.RenderConstants.*;

/**
 * @author DaPorkchop_
 */
@Getter
public abstract class IndexedMultidrawMultipassRenderStrategy<POS extends IFarPos, T extends IFarTile> extends IndexedMultidrawRenderStrategy<POS, T> implements IMultipassRenderStrategy<POS, T> {
    protected final IDrawCommandBuffer[][] passes;

    public IndexedMultidrawMultipassRenderStrategy(int vertexSize) {
        super(vertexSize);

        this.passes = IntStream.range(0, RENDER_PASS_COUNT)
                .mapToObj(pass -> IntStream.range(0, 32).mapToObj(level -> this.createCommandBuffer()).toArray(IDrawCommandBuffer[]::new))
                .toArray(IDrawCommandBuffer[][]::new);
    }

    @Override
    protected void doRelease() {
        super.doRelease();

        for (IDrawCommandBuffer[] pass : this.passes) {
            for (IDrawCommandBuffer buffer : pass) {
                buffer.delete();
            }
        }
    }
}
