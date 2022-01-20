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

#include "parallel.h"

threadPool::threadPool(int n) : threadNumber(n), running(false) {
}

threadPool::~threadPool() {
    if (running)
        stop();
}

void threadPool::start() {
    running = true;
    for (int i = 0; i < threadNumber; i++)
        threads.emplace_back(std::thread(&threadPool::work, this));
}

void threadPool::stop() {
    {
        std::unique_lock<std::mutex> lk(mtx);
        running = false;
        condVar.notify_all();
    }
    for (std::thread &t: threads) {
        if (t.joinable())
            t.join();
    }
}

void threadPool::work() {
    while (running) {
        std::function<void()> task;
        {
            std::unique_lock<std::mutex> lk(mtx);
            if (!tasks.empty()) {
                task = tasks.front();
                tasks.pop();
            } else if (running && tasks.empty())
                condVar.wait(lk);
        }
        if (task)
            task();
    }
}

void threadPool::appendTask(const std::function<void()> &task) {
    if (running) {
        std::unique_lock<std::mutex> lk(mtx);
        tasks.push(task);
        condVar.notify_one();
    }
}