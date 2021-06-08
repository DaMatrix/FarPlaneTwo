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

package net.daporkchop.fp2.mode.voxel.client;

import lombok.NonNull;
import net.daporkchop.fp2.client.DrawMode;
import net.daporkchop.fp2.client.gl.commandbuffer.IDrawCommandBuffer;
import net.daporkchop.fp2.client.gl.commandbuffer.xfb.VanillaTransformFeedbackCommandBuffer;
import net.minecraft.util.BlockRenderLayer;

/**
 * @author DaPorkchop_
 */
public class TransformFeedbackVoxelRenderStrategy extends AbstractIndexedMultidrawVoxelRenderStrategy {
    @Override
    public IDrawCommandBuffer createCommandBuffer() {
        return new VanillaTransformFeedbackCommandBuffer(super.createCommandBuffer(), VoxelShaders.BLOCK_SHADER_TRANSFORM_FEEDBACK);
    }

    @Override
    public void render(@NonNull BlockRenderLayer layer, boolean pre) {
        try (DrawMode mode = DrawMode.LEGACY.begin()) {
            if (layer == BlockRenderLayer.SOLID && !pre) {
                this.preRender();
                this.renderSolid(this.passes[0]);
                this.postRender();
            } else if (layer == BlockRenderLayer.CUTOUT_MIPPED && !pre) {
                this.preRender();
                this.renderCutout(this.passes[1]);
                this.postRender();
            } else if (layer == BlockRenderLayer.TRANSLUCENT && pre) {
                this.preRender();
                this.renderTransparent(this.passes[2]);
                this.postRender();
            }
        }
    }
}
