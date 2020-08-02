/*
 * Adapted from The MIT License (MIT)
 *
 * Copyright (c) 2020-2020 DaPorkchop_
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

package net.daporkchop.fp2.util.threading.executor;

import io.netty.util.concurrent.EventExecutor;
import lombok.NonNull;
import net.daporkchop.lib.concurrent.PExecutors;
import net.daporkchop.lib.concurrent.future.DefaultPFuture;
import net.daporkchop.lib.unsafe.PUnsafe;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Executor;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.function.ObjIntConsumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static net.daporkchop.lib.common.util.PValidation.*;

/**
 * @author DaPorkchop_
 */
public class LazyPriorityExecutor<K extends LazyKey<K>> {
    protected final Thread[] threads;

    protected final Comparator<LazyTask<K, ?, ?>> comparator;
    protected final BlockingQueue<LazyTask<K, ?, ?>> queue;

    protected volatile boolean running;

    public LazyPriorityExecutor(int threads, @NonNull ThreadFactory threadFactory) {
        this.comparator = (a, b) -> a.key() != null ? a.key().compareTo(b.key()) : -1;
        this.queue = new PriorityBlockingQueue<>(256, this.comparator);

        this.threads = new Thread[positive(threads, "threads")];

        Runnable worker = new Worker();
        for (int i = 0; i < threads; i++) {
            (this.threads[i] = threadFactory.newThread(worker)).start();
        }
    }

    public void submit(@NonNull LazyTask<K, ?, ?> task) {
        this.queue.add(task);
    }

    public void submit(@NonNull Collection<LazyTask<K, ?, ?>> tasks) {
        this.queue.addAll(tasks);
    }

    public void shutdown() {
        this.running = false;

        //interrupt all workers
        for (Thread t : this.threads) {
            t.interrupt();
        }

        //wait for all workers to shut down
        for (Thread t : this.threads) {
            do {
                try {
                    t.join();
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }
            } while (t.isAlive());
        }
    }

    public EventExecutor maxPriorityExecutor()   {
        class MaxPriorityExecutor implements Executor {
            protected final EventExecutor nettyExecutor = PExecutors.toNettyExecutor(this);

            @Override
            public void execute(@NonNull Runnable task) {
                class RunnableWrapper extends DefaultPFuture<Void> implements LazyTask<K,Void, Void> {
                    protected final Runnable task;

                    public RunnableWrapper(@NonNull Runnable task) {
                        super(MaxPriorityExecutor.this.nettyExecutor);

                        this.task = task;
                    }

                    @Override
                    public K key() {
                        return null;
                    }

                    @Override
                    public Stream<? extends LazyTask<K, ?, Void>> before(@NonNull K key) {
                        return Stream.empty();
                    }

                    @Override
                    public Void run(@NonNull List<Void> params, @NonNull LazyPriorityExecutor<K> executor) {
                        this.task.run();
                        return null;
                    }

                    @Override
                    public RunnableWrapper setSuccess(Void result) {
                        super.setSuccess(result);
                        return this;
                    }

                    @Override
                    public RunnableWrapper setFailure(Throwable cause) {
                        super.setFailure(cause);
                        return this;
                    }

                    @Override
                    public void cancel() {
                        super.cancel(false);
                    }
                }
                LazyPriorityExecutor.this.submit(new RunnableWrapper(task));
            }
        }
        return new MaxPriorityExecutor().nettyExecutor;
    }

    protected class Worker implements Runnable {
        @Override
        public void run() {
            while (LazyPriorityExecutor.this.running) {
                try {
                    this.runSingle();
                } catch (InterruptedException e) {
                    //gracefully exit on interrupt
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

        protected void runSingle() throws InterruptedException {
            this.runLazy(LazyPriorityExecutor.this.queue.take());
        }

        protected <T, R> void runLazy(@NonNull LazyTask<K, T, R> task) throws InterruptedException {
            try {
                List<T> params = this.runBefore(task.before(task.key().raiseTie()).collect(Collectors.toList()));

                R val = task.run(params, LazyPriorityExecutor.this);
                if (val != null)    {
                    task.setSuccess(val);
                }
            } catch (InterruptedException e) {
                throw e; //rethrow
            } catch (Exception e) {
                e.printStackTrace();
                task.setFailure(task.cause());
            }
        }

        protected <T> List<T> runBefore(@NonNull List<LazyTask<K, ?, T>> tasks) throws Exception {
            if (tasks.isEmpty()) {
                return Collections.emptyList();
            }

            List<T> list = new ArrayList<>(tasks.size());
            for (int i = 0, size = tasks.size(); i < size; i++) {
                list.add(null);
                int _i = i; //aaaaaa
                tasks.get(i).thenAccept(v -> list.set(_i, v)); //store in list
            }

            LazyPriorityExecutor.this.queue.addAll(tasks);
            try {
                do {
                    this.runSingle();
                } while (this.areAnyIncomplete(tasks));
            } catch (Exception e) {
                if (!(e instanceof InterruptedException)) {
                    LazyPriorityExecutor.this.queue.removeAll(tasks.stream()
                            .filter(f -> !f.isDone())
                            .peek(LazyTask::cancel)
                            .collect(Collectors.toList()));
                }
                throw e;
            }
            return list;
        }

        protected <V> boolean areAnyIncomplete(@NonNull List<LazyTask<K, ?, V>> tasks) {
            for (int i = 0, size = tasks.size(); i < size; i++) {
                if (!tasks.get(i).isDone()) {
                    return false;
                }
            }
            return true;
        }
    }
}