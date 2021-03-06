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

package net.daporkchop.fp2.mode.api.server.storage;

import lombok.NonNull;
import net.daporkchop.fp2.mode.api.Compressed;
import net.daporkchop.fp2.mode.api.IFarPos;
import net.daporkchop.fp2.util.IReusablePersistent;

import java.io.Closeable;
import java.io.IOException;
import java.util.function.Consumer;

/**
 * Handles reading and writing of far terrain data.
 *
 * @author DaPorkchop_
 */
public interface IFarStorage<POS extends IFarPos, V extends IReusablePersistent> extends Closeable {
    /**
     * Loads the value at the given position.
     *
     * @param pos the position of the value to load
     * @return the loaded value, or {@code null} if it doesn't exist
     */
    Compressed<POS, V> load(@NonNull POS pos);

    /**
     * Stores the given value at the given position, atomically replacing any existing value.
     *
     * @param pos   the position to save the data at
     * @param value the value to save
     */
    void store(@NonNull POS pos, @NonNull Compressed<POS, V> value);

    /**
     * @return the {@link IFarDirtyTracker} used by this storage
     */
    IFarDirtyTracker<POS> dirtyTracker();

    /**
     * Closes this storage.
     * <p>
     * If write operations are queued, this method will block until they are completed.
     */
    @Override
    void close() throws IOException;
}
