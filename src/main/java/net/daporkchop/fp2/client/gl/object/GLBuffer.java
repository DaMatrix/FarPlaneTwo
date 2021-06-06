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

package net.daporkchop.fp2.client.gl.object;

import io.netty.buffer.ByteBuf;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import lombok.Getter;
import lombok.NonNull;
import net.daporkchop.fp2.util.DirectBufferReuse;
import net.daporkchop.lib.unsafe.PUnsafe;
import org.lwjgl.opengl.GL15;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;

import static java.lang.Math.*;
import static net.daporkchop.fp2.client.gl.OpenGL.*;
import static net.daporkchop.lib.common.util.PValidation.*;
import static org.lwjgl.opengl.GL15.*;
import static org.lwjgl.opengl.GL30.*;

/**
 * @author DaPorkchop_
 */
public class GLBuffer extends GLObject implements IGLBuffer {
    protected static final IntSet ACTIVE_TARGETS = new IntOpenHashSet();

    @Getter
    protected long capacity = -1L;

    protected final int usage;
    protected int target = -1;

    public GLBuffer(int usage) {
        super(glGenBuffers());
        this.usage = usage;
    }

    @Override
    public GLBuffer bind(int target) {
        checkState(this.target < 0, "buffer id=%s is already bound!", this.id);
        checkState(ACTIVE_TARGETS.add(target), "target id=%s is already active!", target);
        glBindBuffer(this.target = target, this.id);
        return this;
    }

    @Override
    public void close() {
        checkState(this.target >= 0, "not bound!");
        checkState(ACTIVE_TARGETS.remove(this.target));
        glBindBuffer(this.target, 0);
        this.target = -1;
    }

    /**
     * Sets the capacity of this buffer.
     * <p>
     * If the buffer already has the requested capacity, nothing is changed. Otherwise, the buffer is resized and its contents are now undefined.
     *
     * @param capacity the new capacity
     */
    public void capacity(long capacity) {
        checkState(this.target >= 0, "not bound!");
        glBufferData(this.target, this.capacity = notNegative(capacity, "capacity"), this.usage);
    }

    /**
     * Sets the capacity of this buffer.
     * <p>
     * Unlike {@link #capacity(long)}, this method will retain as much of the original data as possible.
     * <p>
     * If the new capacity is less than the current capacity, the buffer's contents will be truncated. If greater than the current capacity, the
     * data will be extended with undefined contents.
     *
     * @param capacity the new capacity
     */
    public void resize(long capacity) {
        checkState(this.target >= 0, "not bound!");
        long retainedCapacity = min(this.capacity, notNegative(capacity, "capacity"));

        if (retainedCapacity <= 0L) { //previous capacity was unset, so no data needs to be retained
            this.capacity(capacity);
            return;
        } else if (this.capacity == capacity) { //capacity remains unchanged, nothing to do!
            return;
        }

        long buffer = PUnsafe.allocateMemory(retainedCapacity);
        try {
            //download data to main memory
            this.downloadRange(0L, buffer, retainedCapacity);

            //update capacity
            this.capacity(capacity);

            //re-upload retained data
            this.uploadRange(0L, buffer, retainedCapacity);
        } finally {
            PUnsafe.freeMemory(buffer);
        }
    }

    /**
     * Sets the buffer contents.
     *
     * @param addr the base address of the data to upload
     * @param size the size of the data (in bytes)
     */
    public void upload(long addr, long size) {
        checkState(this.target >= 0, "not bound!");
        glBufferData(this.target, DirectBufferReuse.wrapByte(addr, toInt(size, "size")), this.usage);
        this.capacity = size;
    }

    /**
     * Sets the buffer contents.
     *
     * @param data the data to upload
     */
    public void upload(@NonNull ByteBuf data) {
        checkState(this.target >= 0, "not bound!");
        int readableBytes = data.readableBytes();

        if (data.hasMemoryAddress()) { //fast-track: upload whole buffer contents at once
            glBufferData(this.target, data.internalNioBuffer(data.readerIndex(), readableBytes), this.usage);
        } else { //assume the buffer is a composite (we don't care about heap buffers lol)
            glBufferData(this.target, readableBytes, this.usage);
            ByteBuffer[] nioBuffers = data.nioBuffers(data.readerIndex(), readableBytes);
            long off = 0L;
            for (ByteBuffer nioBuffer : nioBuffers) {
                int size = nioBuffer.remaining();
                glBufferSubData(this.target, off, nioBuffer);
                off += size;
            }
        }
        this.capacity = readableBytes;
    }

    /**
     * Sets the buffer contents.
     *
     * @param data the data to upload
     */
    public void upload(@NonNull ByteBuffer data) {
        checkState(this.target >= 0, "not bound!");
        glBufferData(this.target, data, this.usage);
        this.capacity = data.remaining();
    }

    /**
     * Sets the buffer contents.
     *
     * @param data the data to upload
     */
    public void upload(@NonNull IntBuffer data) {
        checkState(this.target >= 0, "not bound!");
        glBufferData(this.target, data, this.usage);
        this.capacity = data.remaining() * (long) INT_SIZE;
    }

    /**
     * Sets the buffer contents.
     *
     * @param data the data to upload
     */
    public void upload(@NonNull FloatBuffer data) {
        checkState(this.target >= 0, "not bound!");
        glBufferData(this.target, data, this.usage);
        this.capacity = data.remaining() * (long) FLOAT_SIZE;
    }

    @Override
    public void uploadRange(long start, long addr, long size) {
        checkState(this.target >= 0, "not bound!");
        checkRangeLen(this.capacity, start, size);
        glBufferSubData(this.target, start, DirectBufferReuse.wrapByte(addr, toInt(size, "size")));
    }

    @Override
    public void uploadRange(long start, @NonNull ByteBuf data) {
        checkState(this.target >= 0, "not bound!");
        int readableBytes = data.readableBytes();
        checkRangeLen(this.capacity, start, readableBytes);

        if (data.hasMemoryAddress()) { //fast-track: upload whole buffer contents at once
            glBufferSubData(this.target, start, data.internalNioBuffer(data.readerIndex(), readableBytes));
        } else { //assume the buffer is a composite (we don't care about heap buffers lol)
            for (ByteBuffer nioBuffer : data.nioBuffers(data.readerIndex(), readableBytes)) {
                int size = nioBuffer.remaining();
                glBufferSubData(this.target, start, nioBuffer);
                start += size;
            }
        }
    }

    @Override
    public void uploadRange(long start, @NonNull ByteBuffer data) {
        checkState(this.target >= 0, "not bound!");
        checkRangeLen(this.capacity, start, data.remaining());
        glBufferSubData(this.target, start, data);
    }

    @Override
    public void uploadRange(long start, @NonNull IntBuffer data) {
        checkState(this.target >= 0, "not bound!");
        checkRangeLen(this.capacity, start, data.remaining() * (long) INT_SIZE);
        glBufferSubData(this.target, start, data);
    }

    @Override
    public void uploadRange(long start, @NonNull FloatBuffer data) {
        checkState(this.target >= 0, "not bound!");
        checkRangeLen(this.capacity, start, data.remaining() * (long) FLOAT_SIZE);
        glBufferSubData(this.target, start, data);
    }

    @Override
    public void downloadRange(long start, long addr, long size) {
        checkState(this.target >= 0, "not bound!");
        checkRangeLen(this.capacity, start, size);
        glGetBufferSubData(this.target, start, DirectBufferReuse.wrapByte(addr, toInt(size, "size")));
    }

    @Override
    protected Runnable delete(int id) {
        return () -> glDeleteBuffers(id);
    }

    /**
     * Binds this buffer to a shader buffer binding slot.
     *
     * @param target the binding type
     * @param index  the binding index
     */
    public void bindBase(int target, int index) {
        glBindBufferBase(target, index, this.id);
    }

    /**
     * Maps the contents of this buffer into client address space.
     *
     * @param access one of {@link GL15#GL_READ_ONLY}, {@link GL15#GL_WRITE_ONLY} or {@link GL15#GL_READ_WRITE}
     * @return the memory address of the mapped buffer contents
     */
    public long map(int access) {
        checkState(this.target >= 0, "not bound!");
        checkState(this.capacity >= 0L, "capacity is unset!");
        return PUnsafe.pork_directBufferAddress(glMapBuffer(this.target, access, this.capacity, DirectBufferReuse._BYTE));
    }

    /**
     * Unmaps the mapped buffer contents.
     */
    public void unmap() {
        checkState(this.target >= 0, "not bound!");
        glUnmapBuffer(this.target);
    }
}
