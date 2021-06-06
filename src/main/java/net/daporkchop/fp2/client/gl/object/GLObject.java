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

import lombok.Getter;
import net.daporkchop.lib.unsafe.PCleaner;
import net.minecraft.client.Minecraft;

import static net.daporkchop.lib.common.util.PValidation.*;

/**
 * Base representation of an OpenGL object.
 *
 * @author DaPorkchop_
 */
@Getter
public abstract class GLObject implements AutoCloseable {
    protected final PCleaner cleaner;
    protected final int id;

    public GLObject(int id) {
        checkArg(id != 0, "id must be non-zero!");
        this.id = id;

        Runnable r = this.delete(id);
        this.cleaner = PCleaner.cleaner(this, () -> Minecraft.getMinecraft().addScheduledTask(r)); //create cleaner to ensure that the object is deleted later
    }

    @Override
    public abstract void close();

    /**
     * Deletes the object immediately (rather than relying on GC to delete it eventually).
     */
    public void delete() {
        this.cleaner.clean();
    }

    /**
     * Gets a {@link Runnable} which will delete the object with the given ID when run.
     * <p>
     * Assume that the given ID is always a valid object ID, and that the function will only ever be called once from the client thread.
     * <p>
     * The returned {@link Runnable} must not hold a reference to this instance, either directly or indirectly as it will prevent the JVM from garbage-collecting
     * it.
     *
     * @param id the ID of the object to delete
     * @return a {@link Runnable} which will delete the object with the given ID when run
     */
    protected abstract Runnable delete(int id);

    @Override
    public final int hashCode() {
        return super.hashCode();
    }

    @Override
    public final boolean equals(Object obj) {
        return super.equals(obj);
    }
}
