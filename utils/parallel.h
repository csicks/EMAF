/***************************************************************************
 *
 * Authors:    Yuxuan Chen
 *
 * Copyright (C) 2021 Pattern Recognition and Bioinformatics Group, Shanghai Jiao Tong University
 *
 * Licensed under the GNU General Public License v3.0 (see LICENSE for details)
 *
 * All comments concerning this program package may be sent to e-mail address 'yxchen11@sjtu.edu.cn'
 ***************************************************************************/

#ifndef EM_PARALLEL_H
#define EM_PARALLEL_H

/** @file
 * this file contains classes and functions for parallel computing
*/

#include <thread>
#include <mutex>
#include <vector>
#include <atomic>
#include <condition_variable>
#include <queue>
#include <future>
#include <functional>
#include <iostream>

/** @brief
 * class for thread pool of fixed number of threads
 * @example
 * @code
 * threadPool thread_pool(5);
 * thread_pool.start();
 * std::future<int> f[500];
 * for (int i = 0; i < 500; i++) {
 *      f[i] = thread_pool.appendTask(fun, i);
 * }
 * for (auto &i: f) {
 *      cout << i.get() << endl;
 * }
 * thread_pool.stop();
 * @endcode
*/
class threadPool {
public:
    int threadNumber; // thread number of thread pool
    std::atomic_bool running; // whether thread pool is running
    std::mutex mtx; // mutex lock
    std::condition_variable condVar; // condition variable
    std::vector<std::thread> threads; // threads in thread pool
    std::queue<std::function<void()>> tasks; // tasks to be carried out in thread pool

public:
    explicit threadPool(int n);

    /// disable default constructor
    threadPool(const threadPool &tp) = delete;

    /// disable default constructor
    threadPool &operator=(const threadPool &tp) = delete;

    ~threadPool();

    /// start/resume thread pool
    void start();

    /// stop thread pool
    void stop();

    /// work function in each thread
    void work();

    /** append task to thread pool
     * @implements for tasks don't require return values
    */
    void appendTask(const std::function<void()> &task);

    /** append task to thread pool
     * @implements for tasks require return values
     * @example @code
     *              std::future<int> f;
     *              f=thread_pool.appendTask(fun, i);
     *              auto r=f.get();
     *          @endcode
    */
    template<typename F, typename... Args>
    auto appendTask(F &&f, Args &&... args) -> std::future<decltype(f(args...))> {
        using RetType = decltype(f(args...));
        auto task = std::make_shared<std::packaged_task<RetType()>>(
                std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        std::future<RetType> future = task->get_future();
        {
            std::lock_guard<std::mutex> lock{mtx};
            tasks.emplace([task]() {
                (*task)();
            });
        }
        condVar.notify_one();
        return future;
    }

};

#endif //EM_PARALLEL_H
