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

import io.netty.buffer.ByteBuf;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import lombok.NonNull;
import lombok.experimental.UtilityClass;
import net.daporkchop.fp2.client.TexUVs;
import net.daporkchop.fp2.client.gl.object.IGLBuffer;
import net.daporkchop.fp2.client.gl.object.VertexArrayObject;
import net.daporkchop.fp2.client.gl.type.Int2_10_10_10_Rev;
import net.daporkchop.fp2.client.gl.vertex.IVertexAttribute;
import net.daporkchop.fp2.client.gl.vertex.VertexAttributeInterpretation;
import net.daporkchop.fp2.client.gl.vertex.VertexAttributeType;
import net.daporkchop.fp2.client.gl.vertex.VertexFormat;
import net.daporkchop.fp2.compat.vanilla.FastRegistry;
import net.daporkchop.fp2.mode.common.client.BakeOutput;
import net.daporkchop.fp2.mode.voxel.VoxelData;
import net.daporkchop.fp2.mode.voxel.VoxelPos;
import net.daporkchop.fp2.mode.voxel.VoxelTile;
import net.daporkchop.fp2.util.SingleBiomeBlockAccess;
import net.daporkchop.fp2.util.datastructure.PointOctree3I;
import net.daporkchop.lib.common.util.PArrays;
import net.minecraft.block.state.IBlockState;
import net.minecraft.init.Biomes;
import net.minecraft.util.math.BlockPos;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.stream.IntStream;
import java.util.stream.LongStream;
import java.util.stream.Stream;

import static java.lang.Math.*;
import static net.daporkchop.fp2.client.ClientConstants.*;
import static net.daporkchop.fp2.client.gl.GLCompatibilityHelper.*;
import static net.daporkchop.fp2.client.gl.OpenGL.*;
import static net.daporkchop.fp2.mode.voxel.VoxelConstants.*;
import static net.daporkchop.fp2.util.BlockType.*;
import static net.daporkchop.fp2.util.Constants.*;
import static net.daporkchop.fp2.util.math.MathUtil.*;
import static net.daporkchop.lib.common.util.PorkUtil.*;

/**
 * Shared code for baking voxel geometry.
 *
 * @author DaPorkchop_
 */
@UtilityClass
public class VoxelBake {
    protected static final IVertexAttribute.Int1 ATTRIB_STATE = IVertexAttribute.Int1.builder()
            .alignAndPadTo(EFFECTIVE_VERTEX_ATTRIBUTE_ALIGNMENT)
            .type(VertexAttributeType.UNSIGNED_INT)
            .interpretation(VertexAttributeInterpretation.INTEGER)
            .build();

    protected static final IVertexAttribute.Int2 ATTRIB_LIGHT = IVertexAttribute.Int2.builder(ATTRIB_STATE)
            .alignAndPadTo(EFFECTIVE_VERTEX_ATTRIBUTE_ALIGNMENT)
            .type(VertexAttributeType.UNSIGNED_BYTE)
            .interpretation(VertexAttributeInterpretation.NORMALIZED_FLOAT)
            .build();

    protected static final IVertexAttribute.Int3 ATTRIB_COLOR = IVertexAttribute.Int3.builder(ATTRIB_LIGHT)
            .alignAndPadTo(EFFECTIVE_VERTEX_ATTRIBUTE_ALIGNMENT)
            .type(VertexAttributeType.UNSIGNED_BYTE)
            .interpretation(VertexAttributeInterpretation.NORMALIZED_FLOAT)
            .build();

    protected static final IVertexAttribute.Int4 ATTRIB_POS = IVertexAttribute.Int4.builder(ATTRIB_COLOR)
            .alignAndPadTo(EFFECTIVE_VERTEX_ATTRIBUTE_ALIGNMENT)
            .type(WORKAROUND_AMD_INT_2_10_10_10_REV ? VertexAttributeType.SHORT : VertexAttributeType.INT_2_10_10_10_REV)
            .interpretation(VertexAttributeInterpretation.FLOAT)
            .build();

    protected static final VertexFormat VERTEX_FORMAT = new VertexFormat(ATTRIB_POS, max(EFFECTIVE_VERTEX_ATTRIBUTE_ALIGNMENT, INT_SIZE));

    public void vertexAttributes(@NonNull IGLBuffer buffer, @NonNull VertexArrayObject vao) {
        FP2_LOG.info("voxel vertex size: {} bytes", VERTEX_FORMAT.size());
        VERTEX_FORMAT.configureVAO(vao, buffer);
    }

    public void bakeForShaderDraw(@NonNull VoxelPos dstPos, @NonNull VoxelTile[] srcs, @NonNull BakeOutput output, @NonNull ByteBuf verts, @NonNull ByteBuf[] indices) {
        if (srcs[0] == null) {
            return;
        }

        //step 1: write vertices for the source tile, and assign indices
        writeVertices(dstPos, srcs[0], verts);

        //step 2: write indices to actually connect the vertices and build the mesh
        writeIndices(srcs[0], indices);

        //force parent tile to be rendered over this one
        output.forceRenderParent = true;
    }

    protected void writeVertices(@NonNull VoxelPos tilePos, @NonNull VoxelTile tile, @NonNull ByteBuf vertices) {
        SingleBiomeBlockAccess biomeAccess = new SingleBiomeBlockAccess();
        VoxelData data = new VoxelData();
        BlockPos.MutableBlockPos blockPos = new BlockPos.MutableBlockPos();

        final int level = tilePos.level();
        final int blockX = tilePos.blockX();
        final int blockY = tilePos.blockY();
        final int blockZ = tilePos.blockZ();

        for (int i = 0, lim = tile.vertexCount(); i < lim; i++) {
            tile.getVertex(i, data);

            biomeAccess.biome(FastRegistry.getBiome(data.biome, Biomes.PLAINS));
            blockPos.setPos(blockX + (data.x << level >> POS_FRACT_SHIFT), blockY + (data.y << level >> POS_FRACT_SHIFT), blockZ + (data.z << level >> POS_FRACT_SHIFT));

            int vertexBase = VERTEX_FORMAT.appendVertex(vertices);

            IBlockState state = FastRegistry.getBlockState(data.state);
            ATTRIB_STATE.set(vertices, vertexBase, TexUVs.STATEID_TO_INDEXID.get(state));

            int blockLight = data.light & 0xF;
            int skyLight = data.light >> 4;
            ATTRIB_LIGHT.set(vertices, vertexBase, blockLight | (blockLight << 4), skyLight | (skyLight << 4));
            ATTRIB_COLOR.setRGB(vertices, vertexBase, mc.getBlockColors().colorMultiplier(state, biomeAccess, blockPos, 0));

            ATTRIB_POS.setInt2_10_10_10_rev(vertices, vertexBase, Int2_10_10_10_Rev.packXYZ(data.x, data.y, data.z));
        }
    }

    protected void writeIndices(@NonNull VoxelTile tile, @NonNull ByteBuf[] indices) {
        VoxelData data = new VoxelData();

        for (int i = 0, lim = tile.indexCount(); i < lim; i += 3) {
            //the mesh always consists of triangles. we can iterate through one triangle at a time and examine the provoking vertex
            // in order to determine which render layer it should be put on.

            int v0 = tile.getIndex(i + 0);
            int v1 = tile.getIndex(i + 1);
            int v2 = tile.getIndex(i + 2);

            tile.getVertex(v2, data);
            IBlockState state = FastRegistry.getBlockState(data.state);
            indices[renderType(state)].writeShortLE(v0).writeShortLE(v1).writeShortLE(v2).writeShortLE(v2); //need to translate this to a quad
        }
    }
}
