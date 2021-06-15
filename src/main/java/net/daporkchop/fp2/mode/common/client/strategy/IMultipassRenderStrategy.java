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

import lombok.NonNull;
import net.daporkchop.fp2.client.gl.commandbuffer.IDrawCommandBuffer;
import net.daporkchop.fp2.config.FP2Config;
import net.daporkchop.fp2.mode.api.IFarPos;
import net.daporkchop.fp2.mode.api.IFarTile;
import net.daporkchop.fp2.mode.common.client.IFarRenderStrategy;
import net.daporkchop.lib.unsafe.PUnsafe;
import net.minecraft.client.renderer.GlStateManager;
import net.minecraft.client.renderer.texture.TextureMap;

import static net.daporkchop.fp2.client.ClientConstants.*;
import static net.daporkchop.fp2.client.gl.OpenGL.*;
import static net.daporkchop.fp2.debug.FP2Debug.*;
import static net.daporkchop.fp2.mode.common.client.RenderConstants.*;
import static net.daporkchop.lib.common.util.PValidation.*;
import static org.lwjgl.opengl.GL11.*;

/**
 * @author DaPorkchop_
 */
public interface IMultipassRenderStrategy<POS extends IFarPos, T extends IFarTile> extends IFarRenderStrategy<POS, T> {
    /**
     * @return an array containing the {@link IDrawCommandBuffer} used for each render pass
     */
    IDrawCommandBuffer[][] passes();

    @Override
    default void prepareRender(long tilev, int tilec) {
        IDrawCommandBuffer[][] passes = this.passes();
        checkArg(passes.length == RENDER_PASS_COUNT, "invalid number of render passes: %d (expected %d)", passes.length, RENDER_PASS_COUNT);

        for (IDrawCommandBuffer[] pass : passes) { //begin all passes
            for (IDrawCommandBuffer buffer : pass) {
                buffer.begin();
            }
        }
        try {
            for (long end = tilev + (long) LONG_SIZE * tilec; tilev != end; tilev += LONG_SIZE) { //draw each tile individually
                this.drawTile(passes, PUnsafe.getLong(tilev));
            }
        } finally {
            for (IDrawCommandBuffer[] pass : passes) { //finish all passes
                for (IDrawCommandBuffer buffer : pass) {
                    buffer.close();
                }
            }
        }
    }

    void drawTile(@NonNull IDrawCommandBuffer[][] passes, long tile);

    default void preRender() {
        if (FP2_DEBUG && FP2Config.debug.disableBackfaceCull) {
            GlStateManager.disableCull();
        }

        glEnable(GL_STENCIL_TEST);

        glStencilMask(0xFF);
        glClearStencil(0x7F);
        GlStateManager.clear(GL_STENCIL_BUFFER_BIT);
    }

    default void postRender() {
        glDisable(GL_STENCIL_TEST);

        if (FP2_DEBUG && FP2Config.debug.disableBackfaceCull) {
            GlStateManager.enableCull();
        }
    }

    default void renderSolid(@NonNull IDrawCommandBuffer[] draw) {
        GlStateManager.disableAlpha();

        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
        for (int level = 0; level < draw.length; level++) {
            glStencilFunc(GL_LEQUAL, level, 0x7F);
            draw[level].draw();
        }

        GlStateManager.enableAlpha();
    }

    default void renderCutout(@NonNull IDrawCommandBuffer[] draw) {
        mc.getTextureManager().getTexture(TextureMap.LOCATION_BLOCKS_TEXTURE).setBlurMipmap(false, mc.gameSettings.mipmapLevels > 0);

        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
        for (int level = 0; level < draw.length; level++) {
            glStencilFunc(GL_LEQUAL, level, 0x7F);
            draw[level].draw();
        }

        mc.getTextureManager().getTexture(TextureMap.LOCATION_BLOCKS_TEXTURE).restoreLastBlurMipmap();
    }

    default void renderTransparent(@NonNull IDrawCommandBuffer[] draw) {
        this.renderTransparentStencilPass(draw);
        this.renderTransparentFragmentPass(draw);
    }

    default void renderTransparentStencilPass(@NonNull IDrawCommandBuffer[] draw) {
        GlStateManager.colorMask(false, false, false, false);

        GlStateManager.depthMask(false);

        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
        for (int level = 0; level < draw.length; level++) {
            glStencilFunc(GL_GEQUAL, 0x80 | (draw.length - level), 0xFF);
            draw[level].draw();
        }

        GlStateManager.depthMask(true);

        GlStateManager.colorMask(true, true, true, true);
    }

    default void renderTransparentFragmentPass(@NonNull IDrawCommandBuffer[] draw) {
        GlStateManager.enableBlend();
        GlStateManager.tryBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);
        GlStateManager.alphaFunc(GL_GREATER, 0.1f);

        glStencilMask(0);

        glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
        for (int level = 0; level < draw.length; level++) {
            glStencilFunc(GL_EQUAL, 0x80 | (draw.length - level), 0xFF);
            draw[level].draw();
        }

        GlStateManager.disableBlend();
    }
}
