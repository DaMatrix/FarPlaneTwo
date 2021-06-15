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

    protected static final IVertexAttribute.Int3 ATTRIB_POS = IVertexAttribute.Int3.builder(ATTRIB_COLOR)
            .alignAndPadTo(EFFECTIVE_VERTEX_ATTRIBUTE_ALIGNMENT)
            .type(WORKAROUND_AMD_INT_2_10_10_10_REV ? VertexAttributeType.SHORT : VertexAttributeType.INT_2_10_10_10_REV)
            .interpretation(VertexAttributeInterpretation.FLOAT)
            .build();

    protected static final VertexFormat VERTEX_FORMAT = new VertexFormat(ATTRIB_POS, max(EFFECTIVE_VERTEX_ATTRIBUTE_ALIGNMENT, INT_SIZE));

    public void vertexAttributes(@NonNull IGLBuffer buffer, @NonNull VertexArrayObject vao) {
        FP2_LOG.info("voxel vertex size: {} bytes", VERTEX_FORMAT.size());
        VERTEX_FORMAT.configureVAO(vao, buffer);
    }

    protected static int vertexMapIndex(int dx, int dy, int dz, int i, int edge) {
        int j = CONNECTION_INDICES[i];
        int ddx = dx + ((j >> 2) & 1);
        int ddy = dy + ((j >> 1) & 1);
        int ddz = dz + (j & 1);

        return ((ddx * T_VERTS + ddy) * T_VERTS + ddz) * EDGE_COUNT + edge;
    }

    public void bakeForShaderDraw(@NonNull VoxelPos dstPos, @NonNull VoxelTile[] srcs, @NonNull BakeOutput output, @NonNull ByteBuf verts, @NonNull ByteBuf[] indices) {
        if (srcs[0] == null) {
            return;
        }

        final int level = dstPos.level();
        final int blockX = dstPos.blockX();
        final int blockY = dstPos.blockY();
        final int blockZ = dstPos.blockZ();

        /*//step 1: write vertices for all source tiles, and assign indices
        writeVertices(srcs, blockX, blockY, blockZ, level, map, verts);

        //step 2: write vertices for all source tiles, and assign indices
        writeVertices(srcs, blockX, blockY, blockZ, level, lowOctree, highOctree, map, verts, output);

        //step 3: write indices to actually connect the vertices and build the mesh
        writeIndices(srcs[0], map, indices, lowOctree);*/

        writeVertices(dstPos, srcs[0], srcs, verts, buildLowPointOctrees(srcs), buildHighPointOctree(srcs, dstPos));
        writeIndices(srcs[0], indices);
    }

    protected void writeVertices(@NonNull VoxelPos tilePos, @NonNull VoxelTile tile, @NonNull VoxelTile[] srcs, @NonNull ByteBuf vertices, @NonNull PointOctree3I[] octrees, PointOctree3I highOctree) {
        SingleBiomeBlockAccess biomeAccess = new SingleBiomeBlockAccess();
        VoxelData data = new VoxelData();
        BlockPos.MutableBlockPos blockPos = new BlockPos.MutableBlockPos();

        BitSet[] allUnusedSrcVerts = IntStream.range(0, 8).mapToObj(i -> new BitSet()).toArray(BitSet[]::new);
        BitSet[] allUnusedDstVerts = IntStream.range(0, 8).mapToObj(i -> new BitSet()).toArray(BitSet[]::new);

        for (int i = 0, lim = tile.vertexCount(); i < lim; i++) {
            tile.getVertex(i, data);
            allUnusedSrcVerts[data.highEdge].set(i);
        }

        for (int ti = 1, tx = 0; tx <= 1; tx++) {
            for (int ty = 0; ty <= 1; ty++) {
                for (int tz = 1; tz <= 1; tz++, ti++) {
                    VoxelTile t = srcs[ti];
                    if (t == null) {
                        continue;
                    }

                    for (int j = 0, lim = t.vertexCount(); j < lim; j++) {
                        t.getVertex(j, data);
                        if (data.lowEdge == ti) {
                            allUnusedDstVerts[ti].set(j);
                        }
                    }
                }
            }
        }

        List<int[]>[] triangles = buildTriangleIndex(tile);

        for (int ti = 0, tx = 0; tx <= 1; tx++) {
            for (int ty = 0; ty <= 1; ty++) {
                for (int tz = 0; tz <= 1; tz++, ti++) {
                    BitSet unusedSrcVerts = allUnusedSrcVerts[ti];
                    BitSet usedSrcVerts = new BitSet();
                    BitSet unusedDstVerts = allUnusedDstVerts[ti];
                    BitSet usedDstVerts = new BitSet();

                    if (unusedSrcVerts.isEmpty() || unusedDstVerts.isEmpty()) {
                        continue;
                    }

                    int[] srcIndicesReal = new int[unusedSrcVerts.cardinality()];
                    int[] srcPositions = srcIndicesReal.clone();
                    for (int i = -1, j = 0; (i = unusedSrcVerts.nextSetBit(i + 1)) >= 0; j++) {
                        srcIndicesReal[j] = i;

                        tile.getVertex(i, data);
                        srcPositions[j] = Int2_10_10_10_Rev.packXYZ(data.x, data.y, data.z);
                    }

                    int[] dstIndicesReal = new int[unusedDstVerts.cardinality()];
                    int[] dstPositions = dstIndicesReal.clone();
                    for (int i = -1, j = 0; (i = unusedDstVerts.nextSetBit(i + 1)) >= 0; j++) {
                        dstIndicesReal[j] = i;

                        srcs[ti].getVertex(i, data);
                        dstPositions[j] = Int2_10_10_10_Rev.packXYZ(
                                (tx << (T_SHIFT + POS_FRACT_SHIFT)) + data.x,
                                (ty << (T_SHIFT + POS_FRACT_SHIFT)) + data.y,
                                (tz << (T_SHIFT + POS_FRACT_SHIFT)) + data.z);
                    }

                    int[] lengths = new int[srcIndicesReal.length * dstIndicesReal.length];
                    for (int oi = 0, si = 0; si < srcIndicesReal.length; si++) {
                        for (int di = 0; di < dstIndicesReal.length; di++, oi++) {
                            int spos = srcPositions[si];
                            int dpos = dstPositions[di];
                            lengths[oi] =
                                    sq(Int2_10_10_10_Rev.unpackX(spos) - Int2_10_10_10_Rev.unpackX(dpos)) +
                                    sq(Int2_10_10_10_Rev.unpackY(spos) - Int2_10_10_10_Rev.unpackY(dpos)) +
                                    sq(Int2_10_10_10_Rev.unpackZ(spos) - Int2_10_10_10_Rev.unpackZ(dpos));
                        }
                    }

                    Integer[] list = IntStream.range(0, lengths.length).boxed().toArray(Integer[]::new);
                    Arrays.sort(list, Comparator.comparingInt(i -> lengths[i]));

                    IntList[] toIndices = IntStream.range(0, srcIndicesReal.length).mapToObj(i -> new IntArrayList()).toArray(IntList[]::new);
                    for (Integer pi : list) {
                        int si = pi / dstIndicesReal.length;
                        int di = pi % dstIndicesReal.length;

                        if (!usedSrcVerts.get(si) || !usedDstVerts.get(di)) {
                            usedSrcVerts.set(si);
                            usedDstVerts.set(di);
                            toIndices[si].add(di);
                        }
                    }

                    for (int si = 0; si < srcIndicesReal.length; si++) {
                        int realSrcIndex = srcIndicesReal[si];
                        tile.getVertex(realSrcIndex, data);

                        int dstPos = dstPositions[toIndices[si].getInt(0)];
                        data.x = Int2_10_10_10_Rev.unpackX(dstPos);
                        data.y = Int2_10_10_10_Rev.unpackY(dstPos);
                        data.z = Int2_10_10_10_Rev.unpackZ(dstPos);
                        tile.setVertex(realSrcIndex, data);

                        for (int i = 1, lim = toIndices[si].size(); i < lim; i++) {
                            dstPos = dstPositions[toIndices[si].getInt(i)];
                            data.x = Int2_10_10_10_Rev.unpackX(dstPos);
                            data.y = Int2_10_10_10_Rev.unpackY(dstPos);
                            data.z = Int2_10_10_10_Rev.unpackZ(dstPos);
                            int newVertexIndex = tile.appendVertex(data);

                            for (int[] triangle : triangles[realSrcIndex]) {
                                tile.appendTriangle(
                                        triangle[0] == realSrcIndex ? newVertexIndex : triangle[0],
                                        triangle[1] == realSrcIndex ? newVertexIndex : triangle[1],
                                        triangle[2] == realSrcIndex ? newVertexIndex : triangle[2]);
                            }
                        }
                    }

                    /*if (unusedSrcVerts.cardinality() < unusedDstVerts.cardinality()) {
                        continue;
                    }

                    int[] tmpArrSrc = new int[unusedSrcVerts.cardinality()];
                    Int2IntMap srcPosToIndex = new Int2IntOpenHashMap(tmpArrSrc.length);
                    for (int i = -1, j = 0; (i = unusedSrcVerts.nextSetBit(i + 1)) >= 0; j++) {
                        tile.getVertex(i, data);
                        tmpArrSrc[j] = Int2_10_10_10_Rev.packXYZ(data.x, data.y, data.z);
                        srcPosToIndex.put(tmpArrSrc[j], i);
                    }
                    PointOctree3I srcOctree = new PointOctree3I(tmpArrSrc);

                    int[] tmpArrDst = new int[unusedDstVerts.cardinality()];
                    Int2IntMap dstPosToIndex = new Int2IntOpenHashMap(tmpArrDst.length);
                    for (int i = -1, j = 0; (i = unusedDstVerts.nextSetBit(i + 1)) >= 0; j++) {
                        srcs[ti].getVertex(i, data);
                        tmpArrDst[j] = Int2_10_10_10_Rev.packXYZ(
                                (tx << (T_SHIFT + POS_FRACT_SHIFT)) + data.x,
                                (ty << (T_SHIFT + POS_FRACT_SHIFT)) + data.y,
                                (tz << (T_SHIFT + POS_FRACT_SHIFT)) + data.z);
                        dstPosToIndex.put(tmpArrDst[j], i);
                    }
                    PointOctree3I dstOctree = new PointOctree3I(tmpArrDst);

                    for (int i = -1, j = 0; !unusedSrcVerts.isEmpty() && (i = unusedDstVerts.nextSetBit(i + 1)) >= 0; j++) {
                        int dstPos = tmpArrDst[j];
                        int srcPos = srcOctree.nearestNeighborMatching(
                                Int2_10_10_10_Rev.unpackX(dstPos), Int2_10_10_10_Rev.unpackY(dstPos), Int2_10_10_10_Rev.unpackZ(dstPos),
                                (point, x, y, z) -> !usedSrcVerts.get(srcPosToIndex.get(point)));

                        int srcIndex = srcPosToIndex.get(srcPos);
                        usedSrcVerts.set(srcIndex);

                        tile.getVertex(srcIndex, data);
                        data.x = Int2_10_10_10_Rev.unpackX(dstPos);
                        data.y = Int2_10_10_10_Rev.unpackY(dstPos);
                        data.z = Int2_10_10_10_Rev.unpackZ(dstPos);
                        tile.setVertex(srcIndex, data);
                    }*/
                }
            }
        }

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

            ATTRIB_POS.set(vertices, vertexBase, data.x, data.y, data.z);
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
            indices[renderType(state)].writeShortLE(v0).writeShortLE(v1).writeShortLE(v2);
        }
    }

    protected PointOctree3I[] buildLowPointOctrees(VoxelTile[] srcs) {
        PointOctree3I[] out = new PointOctree3I[8];

        final VoxelData data = new VoxelData();
        final IntList lowPoints = new IntArrayList();

        for (int i = 0, tx = 0; tx <= 1; tx++) {
            for (int ty = 0; ty <= 1; ty++) {
                for (int tz = 0; tz <= 1; tz++, i++) {
                    VoxelTile tile = srcs[i];
                    if (tile == null || i == 0) {
                        continue;
                    }

                    for (int j = 0, lim = tile.vertexCount(); j < lim; j++) {
                        tile.getVertex(j, data);

                        int px = (tx << (T_SHIFT + POS_FRACT_SHIFT)) + data.x;
                        int py = (ty << (T_SHIFT + POS_FRACT_SHIFT)) + data.y;
                        int pz = (tz << (T_SHIFT + POS_FRACT_SHIFT)) + data.z;

                        if (px >= Int2_10_10_10_Rev.MIN_XYZ_VALUE && px <= Int2_10_10_10_Rev.MAX_XYZ_VALUE
                            && py >= Int2_10_10_10_Rev.MIN_XYZ_VALUE && py <= Int2_10_10_10_Rev.MAX_XYZ_VALUE
                            && pz >= Int2_10_10_10_Rev.MIN_XYZ_VALUE && pz <= Int2_10_10_10_Rev.MAX_XYZ_VALUE) { //this will only discard a very small minority of vertices
                            lowPoints.add(Int2_10_10_10_Rev.packXYZ(px, py, pz));
                        }
                    }

                    out[i] = lowPoints.isEmpty() ? null : new PointOctree3I(lowPoints.toIntArray());
                    lowPoints.clear();
                }
            }
        }

        return out;
    }

    protected PointOctree3I buildHighPointOctree(VoxelTile[] srcs, VoxelPos pos) {
        final VoxelData data = new VoxelData();
        final IntList highPoints = new IntArrayList();

        int offX = -(pos.x() & 1) << (T_SHIFT + POS_FRACT_SHIFT);
        int offY = -(pos.y() & 1) << (T_SHIFT + POS_FRACT_SHIFT);
        int offZ = -(pos.z() & 1) << (T_SHIFT + POS_FRACT_SHIFT);

        for (int i = 8, tx = BAKE_HIGH_RADIUS_MIN; tx <= BAKE_HIGH_RADIUS_MAX; tx++) {
            for (int ty = BAKE_HIGH_RADIUS_MIN; ty <= BAKE_HIGH_RADIUS_MAX; ty++) {
                for (int tz = BAKE_HIGH_RADIUS_MIN; tz <= BAKE_HIGH_RADIUS_MAX; tz++, i++) {
                    VoxelTile tile = srcs[i];
                    if (tile == null) {
                        continue;
                    }

                    for (int j = 0; j < tile.vertexCount(); j++) {
                        tile.getVertex(j, data);

                        int px = (tx << (T_SHIFT + POS_FRACT_SHIFT)) + (data.x << 1) + offX;
                        int py = (ty << (T_SHIFT + POS_FRACT_SHIFT)) + (data.y << 1) + offY;
                        int pz = (tz << (T_SHIFT + POS_FRACT_SHIFT)) + (data.z << 1) + offZ;

                        if (px >= Int2_10_10_10_Rev.MIN_XYZ_VALUE && px <= Int2_10_10_10_Rev.MAX_XYZ_VALUE
                            && py >= Int2_10_10_10_Rev.MIN_XYZ_VALUE && py <= Int2_10_10_10_Rev.MAX_XYZ_VALUE
                            && pz >= Int2_10_10_10_Rev.MIN_XYZ_VALUE && pz <= Int2_10_10_10_Rev.MAX_XYZ_VALUE) { //this will only discard a very small minority of vertices
                            highPoints.add(Int2_10_10_10_Rev.packXYZ(px, py, pz));
                        }
                    }
                }
            }
        }

        return highPoints.isEmpty() ? null : new PointOctree3I(highPoints.toIntArray());
    }

    protected List<int[]>[] buildTriangleIndex(@NonNull VoxelTile tile) {
        List<int[]>[] lists = uncheckedCast(IntStream.range(0, tile.vertexCount()).mapToObj(i -> new ArrayList<>()).toArray(List[]::new));

        for (int i = 0, lim = tile.triangleCount(); i < lim; i++) {
            int[] triangle = new int[3];
            tile.getTriangle(i, triangle);

            for (int v : triangle) {
                lists[v].add(triangle);
            }
        }

        return lists;
    }
}
