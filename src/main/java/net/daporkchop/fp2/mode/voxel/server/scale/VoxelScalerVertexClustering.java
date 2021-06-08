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

package net.daporkchop.fp2.mode.voxel.server.scale;

import it.unimi.dsi.fastutil.Hash;
import it.unimi.dsi.fastutil.objects.ObjectOpenCustomHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import lombok.NonNull;
import net.daporkchop.fp2.mode.api.server.gen.IFarScaler;
import net.daporkchop.fp2.mode.voxel.VoxelData;
import net.daporkchop.fp2.mode.voxel.VoxelPos;
import net.daporkchop.fp2.mode.voxel.VoxelTile;
import net.daporkchop.fp2.util.datastructure.WritableBVH;
import net.minecraft.util.math.AxisAlignedBB;
import net.minecraft.util.math.Vec3d;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

import static java.lang.Math.*;
import static net.daporkchop.fp2.mode.voxel.VoxelConstants.*;
import static net.daporkchop.fp2.util.Constants.*;
import static net.daporkchop.fp2.util.math.MathUtil.*;
import static net.daporkchop.lib.common.util.PValidation.*;

/**
 * @author DaPorkchop_
 * @see <a href="https://www.researchgate.net/publication/220792145_Model_Simplification_Using_Vertex-Clustering">this paper</a>
 */
public class VoxelScalerVertexClustering implements IFarScaler<VoxelPos, VoxelTile> {
    protected static final int SRC_MIN = 0;
    protected static final int SRC_MAX = T_VOXELS << 1;
    protected static final int SRC_SIZE = SRC_MAX - SRC_MIN;

    protected static final int SRC_TILE_MIN = SRC_MIN >> T_SHIFT;
    protected static final int SRC_TILE_MAX = ((SRC_MAX - 1) >> T_SHIFT) + 1;
    protected static final int SRC_TILE_SIZE = SRC_TILE_MAX - SRC_TILE_MIN;

    protected static final int DST_MIN = -1;
    protected static final int DST_MAX = T_VOXELS + 1;
    protected static final int DST_SIZE = DST_MAX - DST_MIN;

    protected static boolean outOfBounds(VoxelData data) {
        final int min = SRC_MIN << POS_FRACT_SHIFT >> 1;
        final int max = SRC_MAX << POS_FRACT_SHIFT >> 1;
        return ((SRC_MAX | SRC_MAX) & T_MASK) != 0 && (data.x < min || data.x > max || data.y < min || data.y > max || data.z < min || data.z > max);
    }

    protected static BoundedVoxelData findBoundedVoxelData(@NonNull WritableBVH<BoundedVoxelData> bvh, @NonNull VoxelData data) {
        for (BoundedVoxelData boundedData : bvh.intersecting(new Vec3d(data.x, data.y, data.z))) {
            if (((boundedData.data.lowEdge << EDGE_COUNT) | boundedData.data.highEdge) == ((data.lowEdge << EDGE_COUNT) | data.highEdge)) {
                return boundedData;
            }
        }
        throw new IllegalStateException();
    }

    @Override
    public Stream<VoxelPos> outputs(@NonNull VoxelPos srcPos) {
        return Stream.of(srcPos.up()); //TODO: fix this
    }

    @Override
    public Stream<VoxelPos> inputs(@NonNull VoxelPos dstPos) {
        checkArg(dstPos.level() > 0, "cannot generate inputs for level 0!");

        int x = dstPos.x() << 1;
        int y = dstPos.y() << 1;
        int z = dstPos.z() << 1;
        int level = dstPos.level() - 1;

        VoxelPos[] positions = new VoxelPos[cb(SRC_TILE_SIZE)];
        for (int i = 0, dx = SRC_TILE_MIN; dx < SRC_TILE_MAX; dx++) {
            for (int dy = SRC_TILE_MIN; dy < SRC_TILE_MAX; dy++) {
                for (int dz = SRC_TILE_MIN; dz < SRC_TILE_MAX; dz++) {
                    positions[i++] = new VoxelPos(level, x + dx, y + dy, z + dz);
                }
            }
        }
        return Arrays.stream(positions);
    }

    @Override
    public long scale(@NonNull VoxelTile[] srcTiles, @NonNull VoxelTile dst) {
        //read all the triangles from all source tiles, compute their weights and add them to the list
        List<ExtraVoxelData> srcVerts = new ArrayList<>();
        for (int ti = 0, tx = SRC_TILE_MIN; tx < SRC_TILE_MAX; tx++) {
            for (int ty = SRC_TILE_MIN; ty < SRC_TILE_MAX; ty++) {
                for (int tz = SRC_TILE_MIN; tz < SRC_TILE_MAX; tz++, ti++) {
                    VoxelTile tile = srcTiles[ti];
                    if (tile == null) {
                        continue;
                    }

                    int mask = (tx == SRC_TILE_MAX - 1 ? 4 : 0) | (ty == SRC_TILE_MAX - 1 ? 2 : 0) | (tz == SRC_TILE_MAX - 1 ? 1 : 0);

                    for (int i = 0, len = tile.indexCount(); i < len; i += 3) {
                        ExtraVoxelData data0 = new ExtraVoxelData();
                        ExtraVoxelData data1 = new ExtraVoxelData();
                        ExtraVoxelData data2 = new ExtraVoxelData();

                        //load vertices
                        tile.getVertex(tile.getIndex(i + 0), data0);
                        tile.getVertex(tile.getIndex(i + 1), data1);
                        tile.getVertex(tile.getIndex(i + 2), data2);

                        data0.x = data0.x + (tx << (POS_FRACT_SHIFT + T_SHIFT)) >> 1;
                        data0.y = data0.y + (ty << (POS_FRACT_SHIFT + T_SHIFT)) >> 1;
                        data0.z = data0.z + (tz << (POS_FRACT_SHIFT + T_SHIFT)) >> 1;
                        data1.x = data1.x + (tx << (POS_FRACT_SHIFT + T_SHIFT)) >> 1;
                        data1.y = data1.y + (ty << (POS_FRACT_SHIFT + T_SHIFT)) >> 1;
                        data1.z = data1.z + (tz << (POS_FRACT_SHIFT + T_SHIFT)) >> 1;
                        data2.x = data2.x + (tx << (POS_FRACT_SHIFT + T_SHIFT)) >> 1;
                        data2.y = data2.y + (ty << (POS_FRACT_SHIFT + T_SHIFT)) >> 1;
                        data2.z = data2.z + (tz << (POS_FRACT_SHIFT + T_SHIFT)) >> 1;

                        data0.lowEdge &= ~mask;
                        data1.lowEdge &= ~mask;
                        data2.lowEdge &= ~mask;
                        data0.highEdge &= mask;
                        data1.highEdge &= mask;
                        data2.highEdge &= mask;

                        if (outOfBounds(data0) && outOfBounds(data1) && outOfBounds(data2)) {
                            continue;
                        }

                        data0.calcWeight(data1, data2);
                        data1.calcWeight(data0, data2);
                        data2.calcWeight(data0, data1);

                        srcVerts.add(data0);
                        srcVerts.add(data1);
                        srcVerts.add(data2);
                    }
                }
            }
        }

        //sort voxels by weight
        ExtraVoxelData[] srcVertsSorted = srcVerts.toArray(new ExtraVoxelData[0]);
        Arrays.sort(srcVertsSorted);

        //build BVH of intersected tiles
        WritableBVH<BoundedVoxelData> bvh = new WritableBVH<>();
        OUTER_LOOP:
        for (ExtraVoxelData data : srcVertsSorted) {
            for (BoundedVoxelData boundedData : bvh.intersecting(new Vec3d(data.x, data.y, data.z))) {
                if (((boundedData.data.lowEdge << EDGE_COUNT) | boundedData.data.highEdge) == ((data.lowEdge << EDGE_COUNT) | data.highEdge)) {
                    if (boundedData.data.state == 0 && data.state != 0) {
                        boundedData.data.state = data.state;
                        dst.setVertex(boundedData.idx, boundedData.data);
                    }
                    continue OUTER_LOOP;
                }
            }
            bvh.add(new BoundedVoxelData(data, dst.appendVertex(data)));
        }

        ObjectSet<int[]> distinctIndices = new ObjectOpenCustomHashSet<>(new Hash.Strategy<int[]>() {
            @Override
            public int hashCode(int[] o) {
                checkArg(o.length == 3);
                return o[0] + o[1] + o[2];
            }

            @Override
            public boolean equals(int[] a, int[] b) {
                if (b == null) {
                    return false;
                }
                checkArg(a.length == 3 && b.length == 3);

                //order-independent equals
                //(this would be really cool to implement in SIMD)
                return (a[0] == b[0] && a[1] == b[1] && a[2] == b[2])
                       || (a[1] == b[0] && a[2] == b[1] && a[0] == b[2])
                       || (a[2] == b[0] && a[0] == b[1] && a[1] == b[2]);
            }
        });

        //find best voxels to round to, and add to set to ensure they aren't duplicated
        for (int i = 0, lim = srcVerts.size(); i < lim; ) {
            int[] indices = {
                    findBoundedVoxelData(bvh, srcVerts.get(i++)).idx,
                    findBoundedVoxelData(bvh, srcVerts.get(i++)).idx,
                    findBoundedVoxelData(bvh, srcVerts.get(i++)).idx
            };

            if (indices[0] != indices[1] && indices[0] != indices[2] && indices[1] != indices[2]) {
                distinctIndices.add(indices);
            }
        }

        //write indices to output tile
        for (int[] triangle : distinctIndices) {
            dst.appendTriangle(triangle[0], triangle[1], triangle[2]);
        }

        return 0L;
    }

    protected static class ExtraVoxelData extends VoxelData implements Comparable<ExtraVoxelData> {
        public double weight;

        public void calcWeight(VoxelData d1, VoxelData d2) {
            double x1 = this.x - d1.x;
            double y1 = this.y - d1.y;
            double z1 = this.z - d1.z;
            double x2 = this.x - d2.x;
            double y2 = this.y - d2.y;
            double z2 = this.z - d2.z;
            this.weight = -cos(acos((x1 * x2 + y1 * y2 + z1 * z1) * fastInvSqrt((x1 * x1 + y1 * y1 + z1 * z1) * (x2 * x2 + y2 * y2 + z2 * z2))) * 0.5d);
        }

        @Override
        public int compareTo(ExtraVoxelData o) {
            return Double.compare(this.weight, o.weight);
        }
    }

    protected static class BoundedVoxelData extends AxisAlignedBB {
        protected final VoxelData data = new VoxelData();
        protected final int idx;

        public BoundedVoxelData(ExtraVoxelData data, int idx) {
            super(data.x - (1 << POS_FRACT_SHIFT) + 1, data.y - (1 << POS_FRACT_SHIFT) + 1, data.z - (1 << POS_FRACT_SHIFT) + 1,
                    data.x + (1 << POS_FRACT_SHIFT) - 1, data.y + (1 << POS_FRACT_SHIFT) - 1, data.z + (1 << POS_FRACT_SHIFT) - 1);

            this.idx = idx;

            this.data.x = data.x;
            this.data.y = data.y;
            this.data.z = data.z;
            this.data.lowEdge = data.lowEdge;
            this.data.highEdge = data.highEdge;
            this.data.state = data.state;
            this.data.light = data.light;
            this.data.biome = data.biome;
        }
    }
}
