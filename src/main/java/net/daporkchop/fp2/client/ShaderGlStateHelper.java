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

package net.daporkchop.fp2.client;

import lombok.NonNull;
import lombok.experimental.UtilityClass;
import net.daporkchop.fp2.client.gl.MatrixHelper;
import net.daporkchop.fp2.client.gl.object.GLBuffer;
import net.daporkchop.fp2.util.DirectBufferReuse;
import net.daporkchop.lib.common.math.PMath;
import net.daporkchop.lib.unsafe.PUnsafe;
import net.minecraft.client.Minecraft;
import net.minecraft.client.renderer.GlStateManager;
import net.minecraft.entity.Entity;

import java.nio.FloatBuffer;

import static net.daporkchop.fp2.client.gl.OpenGL.*;
import static net.daporkchop.fp2.compat.of.OFHelper.*;
import static net.daporkchop.lib.common.math.PMath.*;
import static net.minecraft.util.math.MathHelper.*;
import static org.lwjgl.opengl.GL11.*;
import static org.lwjgl.opengl.GL15.*;
import static org.lwjgl.opengl.GL31.*;

/**
 * @author DaPorkchop_
 */
@UtilityClass
public class ShaderGlStateHelper {
    private final GLBuffer BUFFER = new GLBuffer(GL_STREAM_DRAW);

    private final int OFFSET_CAMERA = 0;
    private final int SIZE_CAMERA = MAT4_SIZE + VEC4_SIZE + IVEC3_SIZE + VEC3_SIZE;

    private final int OFFSET_FOG = OFFSET_CAMERA + SIZE_CAMERA;
    private final int SIZE_FOG = VEC4_SIZE + 4 * FLOAT_SIZE + INT_SIZE;

    private final int TOTAL_SIZE = OFFSET_FOG + SIZE_FOG;

    private final long DATA = PUnsafe.allocateMemory(TOTAL_SIZE);
    private final long ADDR_CAMERA = DATA + OFFSET_CAMERA;
    private final long ADDR_FOG = DATA + OFFSET_FOG;

    public void update(float partialTicks, @NonNull Minecraft mc) {
        { //camera
            long addr = ADDR_CAMERA;

            //mat4 modelviewprojection
            MatrixHelper.getModelViewProjectionMatrix(addr);
            addr += MAT4_SIZE;

            //vec4 anti_flicker_offset
            PUnsafe.putFloat(addr + 0 * FLOAT_SIZE, 0.0f);
            PUnsafe.putFloat(addr + 1 * FLOAT_SIZE, 0.0f);
            PUnsafe.putFloat(addr + 2 * FLOAT_SIZE, ReversedZ.REVERSED ? 0.0001f : -0.0001f);
            PUnsafe.putFloat(addr + 3 * FLOAT_SIZE, 0.0f);
            addr += VEC4_SIZE;

            Entity entity = mc.getRenderViewEntity();
            double x = entity.lastTickPosX + (entity.posX - entity.lastTickPosX) * partialTicks;
            double y = entity.lastTickPosY + (entity.posY - entity.lastTickPosY) * partialTicks;
            double z = entity.lastTickPosZ + (entity.posZ - entity.lastTickPosZ) * partialTicks;

            //ivec3 position_floor
            PUnsafe.putInt(addr + 0 * INT_SIZE, floorI(x));
            PUnsafe.putInt(addr + 1 * INT_SIZE, floorI(y));
            PUnsafe.putInt(addr + 2 * INT_SIZE, floorI(z));
            addr += IVEC3_SIZE;

            //vec3 position_fract
            PUnsafe.putFloat(addr + 0 * FLOAT_SIZE, (float) frac(x));
            PUnsafe.putFloat(addr + 1 * FLOAT_SIZE, (float) frac(y));
            PUnsafe.putFloat(addr + 2 * FLOAT_SIZE, (float) frac(z));
            addr += VEC3_SIZE;
        }

        { //fog
            long addr = ADDR_FOG;

            //vec4 color (already set externally)
            addr += VEC4_SIZE;

            GlStateManager.FogState fogState = GlStateManager.fogState;
            boolean fogEnabled = fogState.fog.currentState
                                 && (!OF || PUnsafe.getInt(mc.gameSettings, OF_FOGTYPE_OFFSET) != OF_OFF);

            //float density
            PUnsafe.putFloat(addr, fogState.density);
            addr += FLOAT_SIZE;

            //float start
            PUnsafe.putFloat(addr, fogState.start);
            addr += FLOAT_SIZE;

            //float end
            PUnsafe.putFloat(addr, fogState.end);
            addr += FLOAT_SIZE;

            //float scale
            PUnsafe.putFloat(addr, 1.0f / (fogState.end - fogState.start));
            addr += FLOAT_SIZE;

            //int mode
            PUnsafe.putInt(addr, fogEnabled ? fogState.mode : 0);
            addr += INT_SIZE;
        }

        try (GLBuffer buffer = BUFFER.bind(GL_UNIFORM_BUFFER)) { //upload
            buffer.upload(DATA, TOTAL_SIZE);
        }
    }

    public void bind() {
        BUFFER.bindBase(GL_UNIFORM_BUFFER, 0);
    }

    public void updateFogColor(@NonNull FloatBuffer buffer) {
        PUnsafe.copyMemory(PUnsafe.pork_directBufferAddress(buffer) + buffer.position() * FLOAT_SIZE, ADDR_FOG, VEC4_SIZE);
    }
}
