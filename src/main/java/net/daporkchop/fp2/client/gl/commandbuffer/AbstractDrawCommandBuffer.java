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

package net.daporkchop.fp2.client.gl.commandbuffer;

import net.daporkchop.fp2.client.gl.object.GLBuffer;
import net.daporkchop.lib.unsafe.PCleaner;
import net.daporkchop.lib.unsafe.PUnsafe;

import java.util.concurrent.atomic.AtomicLong;

import static net.daporkchop.fp2.util.math.MathUtil.*;
import static net.daporkchop.lib.common.util.PValidation.*;
import static org.lwjgl.opengl.GL15.*;

/**
 * @author DaPorkchop_
 */
public abstract class AbstractDrawCommandBuffer implements IDrawCommandBuffer {
    protected final AtomicLong cleanedAddr = new AtomicLong();
    protected final long entrySize;

    protected long addr;
    protected int size;
    protected int capacity;

    protected final GLBuffer buffer = new GLBuffer(GL_STREAM_DRAW);

    public AbstractDrawCommandBuffer(long entrySize) {
        this(entrySize, 1);
    }

    public AbstractDrawCommandBuffer(long entrySize, int initialCapacity) {
        this.entrySize = positive(entrySize, "entrySize");
        this.capacity = positive(initialCapacity, "initialCapacity");
        this.cleanedAddr.set(this.addr = PUnsafe.allocateMemory(this.capacity * this.entrySize));
        PCleaner.cleaner(this, this.cleanedAddr);
    }

    @Override
    public IDrawCommandBuffer begin() {
        this.size = 0;
        return this;
    }

    protected int next() {
        int size = this.size++;
        if (size == this.capacity) {
            this.grow();
        }
        return size;
    }

    protected void grow() {
        this.capacity = this.capacity + asrCeil(this.capacity, 1); // * 1.5
        this.cleanedAddr.set(this.addr = PUnsafe.reallocateMemory(this.addr, this.capacity * this.entrySize));
    }

    @Override
    public void close() {
        if (this.size != 0) {
            this.upload(this.addr, this.size * this.entrySize);
        }
    }

    protected abstract void upload(long addr, long size);

    @Override
    public void draw() {
        if (this.size != 0) {
            this.draw0();
        }
    }

    protected abstract void draw0();

    @Override
    public void delete() {
        this.buffer.delete();
    }
}
