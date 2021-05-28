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

package net.daporkchop.fp2.mode.voxel.server.gen;

import lombok.NonNull;
import net.daporkchop.fp2.mode.common.server.gen.AbstractFarGenerator;
import net.daporkchop.fp2.mode.voxel.VoxelData;
import net.daporkchop.fp2.mode.voxel.VoxelTile;
import net.daporkchop.fp2.util.math.Vector3d;
import net.daporkchop.fp2.util.math.qef.QefSolver;
import net.daporkchop.lib.common.ref.Ref;
import net.daporkchop.lib.common.ref.ThreadRef;
import net.minecraft.world.WorldServer;

import java.util.Arrays;

import static java.lang.Math.*;
import static net.daporkchop.fp2.mode.voxel.VoxelConstants.*;
import static net.daporkchop.fp2.util.BlockType.*;
import static net.daporkchop.fp2.util.Constants.*;
import static net.daporkchop.fp2.util.math.MathUtil.*;
import static net.daporkchop.lib.common.math.PMath.*;
import static net.daporkchop.lib.common.util.PValidation.*;

/**
 * @author DaPorkchop_
 */
public abstract class AbstractVoxelGenerator<PARAM> extends AbstractFarGenerator {
    public static final int DMAP_MIN = 0;
    public static final int DMAP_MAX = T_VERTS + 1;
    public static final int DMAP_SIZE = DMAP_MAX - DMAP_MIN;

    protected static final int DI_ADD_000 = densityIndex(DMAP_MIN + 0, DMAP_MIN + 0, DMAP_MIN + 0);
    protected static final int DI_ADD_001 = densityIndex(DMAP_MIN + 0, DMAP_MIN + 0, DMAP_MIN + 1);
    protected static final int DI_ADD_010 = densityIndex(DMAP_MIN + 0, DMAP_MIN + 1, DMAP_MIN + 0);
    protected static final int DI_ADD_011 = densityIndex(DMAP_MIN + 0, DMAP_MIN + 1, DMAP_MIN + 1);
    protected static final int DI_ADD_100 = densityIndex(DMAP_MIN + 1, DMAP_MIN + 0, DMAP_MIN + 0);
    protected static final int DI_ADD_101 = densityIndex(DMAP_MIN + 1, DMAP_MIN + 0, DMAP_MIN + 1);
    protected static final int DI_ADD_110 = densityIndex(DMAP_MIN + 1, DMAP_MIN + 1, DMAP_MIN + 0);
    protected static final int DI_ADD_111 = densityIndex(DMAP_MIN + 1, DMAP_MIN + 1, DMAP_MIN + 1);

    protected static final int[] DI_ADD = {
            DI_ADD_000, DI_ADD_001,
            DI_ADD_010, DI_ADD_011,
            DI_ADD_100, DI_ADD_101,
            DI_ADD_110, DI_ADD_111
    };

    protected static int densityIndex(int x, int y, int z) {
        return ((x - DMAP_MIN) * DMAP_SIZE + y - DMAP_MIN) * DMAP_SIZE + z - DMAP_MIN;
    }

    protected final Ref<double[][]> densityMapCache = ThreadRef.soft(() -> new double[2][cb(DMAP_SIZE)]);
    protected final Ref<byte[]> typeMapCache = ThreadRef.soft(() -> new byte[cb(DMAP_SIZE)]);
    protected final Ref<VoxelAlignedMeshAssembler> assemblerCache = ThreadRef.soft(VoxelAlignedMeshAssembler::new);

    public AbstractVoxelGenerator(@NonNull WorldServer world) {
        super(world);
    }

    protected byte[] populateTypeMapFromDensityMap(@NonNull double[][] densityMap) {
        byte[] typeMap = this.typeMapCache.get();
        double[] densityMap0 = densityMap[0];
        double[] densityMap1 = densityMap[1];

        //range check here to allow JIT to avoid range checking inside the loop
        checkArg(typeMap.length >= cb(DMAP_SIZE) && densityMap0.length >= cb(DMAP_SIZE) && densityMap1.length >= cb(DMAP_SIZE));

        for (int i = 0; i < cb(DMAP_SIZE); i++) { //set each type flag depending on the value in the corresponding density map layer
            byte type = 0;
            if (densityMap0[i] > 0.0d) {
                type |= BLOCK_TYPE_TRANSPARENT;
            }
            if (densityMap1[i] > 0.0d) {
                type |= BLOCK_TYPE_OPAQUE;
            }
            typeMap[i] = type;
        }

        return typeMap;
    }

    protected void dualContour(int baseX, int baseY, int baseZ, int level, VoxelTile tile, double[][] densityMap, PARAM param) {
        QefSolver qef = new QefSolver();
        VoxelData data = new VoxelData();
        Vector3d vec = new Vector3d();

        //use bit flags to identify voxel types rather than reading from the density map each time to keep innermost loop head tight and cache-friendly
        byte[] tMap = this.populateTypeMapFromDensityMap(densityMap);

        VoxelAlignedMeshAssembler assembler = this.assemblerCache.get().reset();
        int[] tmpStates = new int[EDGE_COUNT];

        for (int dx = 0; dx < T_VERTS; dx++) {
            for (int dy = 0; dy < T_VERTS; dy++) {
                for (int dz = 0; dz < T_VERTS; dz++) {
                    int diBase = densityIndex(dx, dy, dz);

                    //check for intersection data for each corner
                    int corners = 0;
                    for (int i = 0; i < 8; i++) {
                        int di = diBase + DI_ADD[i];
                        corners |= (tMap[di] & 0xFF) << (i << 1);
                    }

                    if (corners == 0 || corners == 0x5555 || corners == 0xAAAA || corners == 0xFFFF) { //if all corners are the same type, this voxel can be safely skipped
                        continue;
                    }

                    Arrays.fill(tmpStates, -1);

                    double totalNx = 0.0d;
                    double totalNy = 0.0d;
                    double totalNz = 0.0d;

                    //populate the QEF with data
                    qef.reset();
                    int edgeCount = 0;
                    int edges = 0;
                    for (int edge = 0; edge < QEF_EDGE_COUNT; edge++) {
                        int c0 = QEF_EDGE_VERTEX_MAP[edge << 1];
                        int c1 = QEF_EDGE_VERTEX_MAP[(edge << 1) | 1];

                        //determine which layer the two corners that make up this edge are on
                        int layer0 = (corners >> (c0 << 1)) & 3;
                        int layer1 = (corners >> (c1 << 1)) & 3;
                        if (layer0 == layer1 //both corners along the current edge are identical, this edge can be skipped
                            || max(layer0, layer1) == 3 && min(layer0, layer1) == 2) { //don't consider edges that transition from solid+liquid to solid, because then we'll be
                            // generating tons of internal mesh for no reason
                            continue;
                        }

                        //the connection will be made on the bottommost layer
                        int layer = min(max(layer0, layer1) - 1, 1);

                        //compute gradients in the current voxel on the chosen layer, which will also serve as the normal vector
                        double densityBase = densityMap[layer][diBase + DI_ADD_000];
                        double nx = densityMap[layer][diBase + DI_ADD_100] - densityBase;
                        double ny = densityMap[layer][diBase + DI_ADD_010] - densityBase;
                        double nz = densityMap[layer][diBase + DI_ADD_001] - densityBase;
                        totalNx += nx;
                        totalNy += ny;
                        totalNz += nz;

                        //get density values at both corners
                        double density0 = densityMap[layer][diBase + DI_ADD[c0]];
                        double density1 = densityMap[layer][diBase + DI_ADD[c1]];

                        //minimize the density function along the edge to find the point where it crosses 0
                        double t = clamp(minimize(density0, density1), 0.0d, 1.0d);
                        double px = lerp((c0 >> 2) & 1, (c1 >> 2) & 1, t);
                        double py = lerp((c0 >> 1) & 1, (c1 >> 1) & 1, t);
                        double pz = lerp(c0 & 1, c1 & 1, t);

                        //add the edge crossing point to the QEF
                        qef.add(px, py, pz, nx, ny, nz);
                        edgeCount++;

                        if ((edge & 3) == 3) { //this is a renderable edge, so we need to set the state and face direction
                            int faceEdge = edge >> 2;
                            if (layer0 < layer1) { //the face is facing towards negative coordinates
                                edges |= EDGE_DIR_NEGATIVE << (faceEdge << 1);
                            } else {
                                edges |= EDGE_DIR_POSITIVE << (faceEdge << 1);
                            }

                            int state = this.getFaceState(baseX + (dx << level), baseY + (dy << level), baseZ + (dz << level), level, nx, ny, nz, density0, density1, faceEdge, layer, param);
                            //data.states[faceEdge] = state;
                            tmpStates[faceEdge] = state;
                        }
                    }

                    //yet another sanity check: a few voxels will make it through the check before the QEF initialization loop (specifically opaque+transparent -> opaque transitions),
                    // so provide the option to break out here without setting the voxel if we can
                    if (edgeCount == 0) {
                        continue;
                    }

                    //for some reason y is backwards, so we invert it as long as it isn't facing both directions
                    if ((((edges >> 2) ^ (edges >> 3)) & 1) != 0) {
                        edges ^= EDGE_DIR_MASK << 2;
                    }

                    //solve QEF and set the tile data
                    qef.solve(vec, 0.1, 1, 0.5);
                    if (vec.x < 0.0d || vec.x > 1.0d
                        || vec.y < 0.0d || vec.y > 1.0d
                        || vec.z < 0.0d || vec.z > 1.0d) { //ensure that all points are within voxel bounds
                        //if not, fall back to the mass point (basically the average position of all edge intersections), which is basically guaranteed to be within the voxel bounds
                        vec.set(qef.massPoint().x, qef.massPoint().y, qef.massPoint().z);
                    }

                    data.x = (dx << POS_FRACT_SHIFT) + clamp(floorI(vec.x * POS_ONE), 0, POS_ONE);
                    data.y = (dy << POS_FRACT_SHIFT) + clamp(floorI(vec.y * POS_ONE), 0, POS_ONE);
                    data.z = (dz << POS_FRACT_SHIFT) + clamp(floorI(vec.z * POS_ONE), 0, POS_ONE);

                    //normalize normal vector
                    double nFactor = 1.0d / sqrt(totalNx * totalNx + totalNy * totalNy + totalNz * totalNz);
                    totalNx *= nFactor;
                    totalNy *= nFactor;
                    totalNz *= nFactor;

                    this.populateVoxelBlockData(baseX + (dx << level), baseY + (dy << level), baseZ + (dz << level), level, totalNx, totalNy, totalNz, data, param);

                    //let assembler write vertex data to tile
                    assembler.setEdgesAndVertices(tile, dx, dy, dz, edges, data, tmpStates);
                }
            }
        }

        assembler.assemble(tile);
    }

    protected abstract int getFaceState(int blockX, int blockY, int blockZ, int level, double nx, double ny, double nz, double density0, double density1, int edge, int layer, PARAM param);

    protected abstract void populateVoxelBlockData(int blockX, int blockY, int blockZ, int level, double nx, double ny, double nz, VoxelData data, PARAM param);
}
