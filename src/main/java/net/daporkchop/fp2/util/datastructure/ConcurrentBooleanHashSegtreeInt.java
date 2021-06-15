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

import lombok.EqualsAndHashCode;
import lombok.NonNull;
import lombok.RequiredArgsConstructor;
import net.daporkchop.fp2.mode.common.client.FarRenderTree;
import net.daporkchop.lib.common.function.throwing.ESupplier;

import java.util.Set;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static net.daporkchop.lib.common.util.PValidation.*;

/**
 * @author DaPorkchop_
 */
//TODO: make this more memory-efficient
public class ConcurrentBooleanHashSegtreeInt {
    protected final Set<Key> set = ConcurrentHashMap.newKeySet();
    protected final int dims;

    protected CompletableFuture<Void> initializationFuture = new CompletableFuture<>();

    public ConcurrentBooleanHashSegtreeInt(int dims, @NonNull ESupplier<Stream<int[]>> coords) {
        this.dims = positive(dims, "dims");

        CompletableFuture.runAsync(() -> {
            try {
                //insert all coordinates into this segtree
                coords.get().parallel().forEach(this::set);

                //complete the future, then set it to null
                this.initializationFuture.complete(null);
                this.initializationFuture = null;
            } catch (Exception e) {
                //complete the future with an exception
                this.initializationFuture.completeExceptionally(e);
            }
        });
    }

    protected Key makeKey(int level, @NonNull int... coords) {
        checkArg(coords.length == this.dims, "expected %d-dimensional array, but found %d!", this.dims, coords.length);

        switch (coords.length) {
            case 2:
                return new Key2(level, coords[0], coords[1]);
            case 3:
                return new Key3(level, coords[0], coords[1], coords[2]);
            default:
                return new KeyAny(level, coords);
        }
    }

    public void set(@NonNull int... coords) {
        if (this.initializationFuture != null && this.initializationFuture.isCompletedExceptionally()) {
            //throw an exception if initialization failed
            this.initializationFuture.join();
        }

        Key key = this.makeKey(0, coords);
        for (int i = 0; i < FarRenderTree.DEPTH && this.set.add(key); i++, key = key.up()) {
        }
    }

    public boolean isSet(int level, @NonNull int... coords) {
        if (this.initializationFuture != null) {
            //join future to wait until this instance is initialized (or throw an exception if initialization failed)
            this.initializationFuture.join();
        }

        return this.set.contains(this.makeKey(level, coords));
    }

    protected static abstract class Key {
        protected abstract Key up();

        @Override
        public abstract int hashCode();

        @Override
        public abstract boolean equals(Object obj);
    }

    @RequiredArgsConstructor
    @EqualsAndHashCode(callSuper = false)
    protected static final class Key2 extends Key {
        protected final int level;
        protected final int x;
        protected final int z;

        @Override
        protected Key up() {
            return new Key2(this.level + 1, this.x >> 1, this.z >> 1);
        }
    }

    @RequiredArgsConstructor
    @EqualsAndHashCode(callSuper = false)
    protected static final class Key3 extends Key {
        protected final int level;
        protected final int x;
        protected final int y;
        protected final int z;

        @Override
        protected Key up() {
            return new Key3(this.level + 1, this.x >> 1, this.y >> 1, this.z >> 1);
        }
    }

    @RequiredArgsConstructor
    @EqualsAndHashCode(callSuper = false)
    protected static final class KeyAny extends Key {
        protected final int level;
        @NonNull
        protected final int[] coords;

        @Override
        public Key up() {
            return new KeyAny(this.level + 1, IntStream.of(this.coords).map(i -> i >> 1).toArray());
        }
    }
}
