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

package net.daporkchop.fp2.util.datastructure;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import lombok.NonNull;
import lombok.RequiredArgsConstructor;
import net.daporkchop.fp2.client.gl.type.Int2_10_10_10_Rev;
import net.daporkchop.fp2.util.math.geometry.Sphere;
import net.daporkchop.fp2.util.math.geometry.Volume;
import net.daporkchop.lib.common.util.PArrays;
import net.minecraft.util.math.AxisAlignedBB;

import java.util.Arrays;

import static java.lang.Math.*;
import static net.daporkchop.fp2.util.math.MathUtil.*;
import static net.daporkchop.lib.common.util.PValidation.*;

/**
 * Immutable, 3-dimensional implementation of an octree over a cloud of points with {@code int} coordinates.
 *
 * @author DaPorkchop_
 */
//TODO: this is a simple implementation, but certainly not particularly fast. let's optimize this to use flattened structures (probably off-heap)
public final class PointOctree3I {
    protected static final int TARGET_POINTS_PER_NODE = 256;

    protected final Node root;

    public PointOctree3I(@NonNull int... points) {
        this.root = points.length > 0 ? new Node(points) : null;
    }

    /**
     * Finds the point in the octree which has the shortest distance to the given point.
     *
     * @param x the point's X coordinate
     * @param y the point's Y coordinate
     * @param z the point's Z coordinate
     * @return the closest point to the given point, or {@code -1} if no points match
     */
    public int nearestNeighbor(int x, int y, int z) {
        Int2_10_10_10_Rev.checkXYZ(x, "x");
        Int2_10_10_10_Rev.checkXYZ(y, "y");
        Int2_10_10_10_Rev.checkXYZ(z, "z");

        if (this.root == null) { //octree is empty! no points match
            return -1;
        }

        NearestNeighborQuery query = new NearestNeighborQuery(x, y, z);
        this.root.nearestNeighbor(query);
        return query.bestNeighbor;
    }

    /**
     * Finds the point in the octree which has the shortest distance to the given point.
     *
     * @param x the point's X coordinate
     * @param y the point's Y coordinate
     * @param z the point's Z coordinate
     * @return the closest point to the given point, or {@code -1} if no points match
     */
    public int nearestNeighborInBounds(int x, int y, int z, int minX, int minY, int minZ, int maxX, int maxY, int maxZ) {
        Int2_10_10_10_Rev.checkXYZ(x, "x");
        Int2_10_10_10_Rev.checkXYZ(y, "y");
        Int2_10_10_10_Rev.checkXYZ(z, "z");

        if (this.root == null) { //octree is empty! no points match
            return -1;
        }

        BoundedNearestNeighborQuery query = new BoundedNearestNeighborQuery(x, y, z, minX, minY, minZ, maxX, maxY, maxZ);
        this.root.nearestNeighbor(query);
        return query.bestNeighbor;
    }

    /**
     * Finds the point in the octree which has the shortest distance to the given point.
     *
     * @param x the point's X coordinate
     * @param y the point's Y coordinate
     * @param z the point's Z coordinate
     * @return the closest point to the given point, or {@code -1} if no points match
     */
    public int nearestNeighborOutsideBounds(int x, int y, int z, int minX, int minY, int minZ, int maxX, int maxY, int maxZ) {
        Int2_10_10_10_Rev.checkXYZ(x, "x");
        Int2_10_10_10_Rev.checkXYZ(y, "y");
        Int2_10_10_10_Rev.checkXYZ(z, "z");

        if (this.root == null) { //octree is empty! no points match
            return -1;
        }

        ExclusiveBoundedNearestNeighborQuery query = new ExclusiveBoundedNearestNeighborQuery(x, y, z, minX, minY, minZ, maxX, maxY, maxZ);
        this.root.nearestNeighbor(query);
        return query.bestNeighbor;
    }

    /**
     * Finds the point in the octree which has the shortest distance to the given point.
     *
     * @param x the point's X coordinate
     * @param y the point's Y coordinate
     * @param z the point's Z coordinate
     * @return the closest point to the given point, or {@code -1} if no points match
     */
    public int nearestNeighborMatching(int x, int y, int z, @NonNull PointPredicate condition) {
        Int2_10_10_10_Rev.checkXYZ(x, "x");
        Int2_10_10_10_Rev.checkXYZ(y, "y");
        Int2_10_10_10_Rev.checkXYZ(z, "z");

        if (this.root == null) { //octree is empty! no points match
            return -1;
        }

        ConditionalNearestNeighborQuery query = new ConditionalNearestNeighborQuery(x, y, z, condition);
        this.root.nearestNeighbor(query);
        return query.bestNeighbor;
    }

    protected static class Node {
        protected final Node[] children;
        protected final int[] points;

        protected final AxisAlignedBB bb;
        protected final int middleX;
        protected final int middleY;
        protected final int middleZ;

        public Node(@NonNull int... points) {
            checkArg(points.length > 0, "at least one point must be given!");

            //compute point cloud bounds and centers
            int minX = Integer.MAX_VALUE;
            int minY = Integer.MAX_VALUE;
            int minZ = Integer.MAX_VALUE;
            int maxX = Integer.MIN_VALUE;
            int maxY = Integer.MIN_VALUE;
            int maxZ = Integer.MIN_VALUE;
            for (int point : points) {
                minX = min(minX, Int2_10_10_10_Rev.unpackX(point));
                maxX = max(maxX, Int2_10_10_10_Rev.unpackX(point));
                minY = min(minY, Int2_10_10_10_Rev.unpackY(point));
                maxY = max(maxY, Int2_10_10_10_Rev.unpackY(point));
                minZ = min(minZ, Int2_10_10_10_Rev.unpackZ(point));
                maxZ = max(maxZ, Int2_10_10_10_Rev.unpackZ(point));
            }
            this.bb = new AxisAlignedBB(minX, minY, minZ, maxX, maxY, maxZ);
            this.middleX = (minX + maxX) >> 1;
            this.middleY = (minY + maxY) >> 1;
            this.middleZ = (minZ + maxZ) >> 1;

            //split points up across all 8 child nodes
            IntList[] buffers = new IntList[8];
            PArrays.fill(buffers, IntArrayList::new);

            for (int point : points) {
                int bdx = (Int2_10_10_10_Rev.unpackX(point) - this.middleX) >>> 31;
                int bdy = (Int2_10_10_10_Rev.unpackY(point) - this.middleY) >>> 31;
                int bdz = (Int2_10_10_10_Rev.unpackZ(point) - this.middleZ) >>> 31;
                buffers[(bdx << 2) | (bdy << 1) | bdz].add(point);
            }

            if (points.length <= TARGET_POINTS_PER_NODE //there's a sufficiently small number of points for this node to be a leaf
                || Arrays.stream(buffers).filter(list -> !list.isEmpty()).count() == 1L) { //all points fit into a single sub-volume?!?
                this.children = null;
                this.points = points;
            } else {
                this.children = new Node[8];
                for (int i = 0; i < 8; i++) {
                    this.children[i] = buffers[i].isEmpty() ? null : new Node(buffers[i].toIntArray());
                }

                this.points = null;
            }
        }

        /**
         * @return whether or not this node is a leaf
         */
        protected boolean isLeaf() {
            return this.children == null;
        }

        protected void nearestNeighbor(@NonNull NearestNeighborQuery query) {
            if (!query.boundingVolume().intersects(this.bb)) { //query doesn't intersect this node, there's nothing to do!
                return;
            }

            if (this.isLeaf()) { //this is a leaf node: do a brute-force search against every point
                for (int point : this.points) {
                    query.update(point, Int2_10_10_10_Rev.unpackX(point), Int2_10_10_10_Rev.unpackY(point), Int2_10_10_10_Rev.unpackZ(point));
                }
                return;
            }

            //step 1: recurse down as far as possible into nodes which are as close as possible to the query point

            //find the child node with the minimum distance
            int bdx = (query.x - this.middleX) >>> 31;
            int bdy = (query.y - this.middleY) >>> 31;
            int bdz = (query.z - this.middleZ) >>> 31;
            Node closestChild = this.children[(bdx << 2) | (bdy << 1) | bdz];

            //recurse into close child
            if (closestChild != null) {
                closestChild.nearestNeighbor(query);
            }

            //step 2: attempt query against all remaining children that still intersect the best query volume

            for (Node child : this.children) {
                if (child != null && child != closestChild && query.boundingVolume().intersects(child.bb)) { //don't query same child twice
                    child.nearestNeighbor(query);
                }
            }
        }
    }

    @RequiredArgsConstructor
    protected static class NearestNeighborQuery {
        protected final int x;
        protected final int y;
        protected final int z;

        protected int bestNeighbor;
        protected int bestNeighborDistSq = Integer.MAX_VALUE;

        public void update(int point, int x, int y, int z) {
            int distanceSq = sq(this.x - x) + sq(this.y - y) + sq(this.z - z);
            if (distanceSq < this.bestNeighborDistSq) {
                this.bestNeighborDistSq = distanceSq;
                this.bestNeighbor = point;
            }
        }

        public Volume boundingVolume() {
            return new Sphere(this.x, this.y, this.z, sqrt(this.bestNeighborDistSq));
        }
    }

    protected static class BoundedNearestNeighborQuery extends NearestNeighborQuery {
        protected final int minX;
        protected final int maxX;
        protected final int minY;
        protected final int maxY;
        protected final int minZ;
        protected final int maxZ;

        public BoundedNearestNeighborQuery(int x, int y, int z, int minX, int minY, int minZ, int maxX, int maxY, int maxZ) {
            super(x, y, z);

            this.minX = minX;
            this.minY = minY;
            this.minZ = minZ;
            this.maxX = maxX;
            this.maxY = maxY;
            this.maxZ = maxZ;
        }

        @Override
        public void update(int point, int x, int y, int z) {
            if (x >= this.minX && x <= this.maxX
                && y >= this.minY && y <= this.maxY
                && z >= this.minZ && z <= this.maxZ) {
                super.update(point, x, y, z);
            }
        }
    }

    protected static class ExclusiveBoundedNearestNeighborQuery extends NearestNeighborQuery {
        protected final int minX;
        protected final int maxX;
        protected final int minY;
        protected final int maxY;
        protected final int minZ;
        protected final int maxZ;

        public ExclusiveBoundedNearestNeighborQuery(int x, int y, int z, int minX, int minY, int minZ, int maxX, int maxY, int maxZ) {
            super(x, y, z);

            this.minX = minX;
            this.minY = minY;
            this.minZ = minZ;
            this.maxX = maxX;
            this.maxY = maxY;
            this.maxZ = maxZ;
        }

        @Override
        public void update(int point, int x, int y, int z) {
            if (!(x >= this.minX && x <= this.maxX
                  && y >= this.minY && y <= this.maxY
                  && z >= this.minZ && z <= this.maxZ)) {
                super.update(point, x, y, z);
            }
        }
    }

    protected static class ConditionalNearestNeighborQuery extends NearestNeighborQuery {
        protected final PointPredicate condition;

        public ConditionalNearestNeighborQuery(int x, int y, int z, @NonNull PointPredicate condition) {
            super(x, y, z);

            this.condition = condition;
        }

        @Override
        public void update(int point, int x, int y, int z) {
            if (this.condition.test(point, x, y, z)) {
                super.update(point, x, y, z);
            }
        }
    }

    @FunctionalInterface
    public interface PointPredicate {
        boolean test(int point, int x, int y, int z);
    }
}
