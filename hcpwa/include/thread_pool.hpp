#include <condition_variable>
#include <cstddef>
#include <functional>
#include <future>
#include <iostream>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

class ThreadPool {
 public:
  explicit ThreadPool(std::size_t num_threads)
      : stop_(false) {
    if (num_threads == 0) {
      throw std::invalid_argument("ThreadPool size must be greater than 0");
    }

    workers_.reserve(num_threads);
    for (std::size_t i = 0; i < num_threads; ++i) {
      workers_.emplace_back([this]() {
        workerLoop();
      });
    }
  }

  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator=(const ThreadPool&) = delete;

  ~ThreadPool() {
    shutdown();
  }

  template <typename F, typename... Args>
  auto enqueue(F&& f, Args&&... args)
      -> std::future<std::invoke_result_t<F, Args...>> {
    using ReturnType = std::invoke_result_t<F, Args...>;

    auto task_ptr = std::make_shared<std::packaged_task<ReturnType()>>(
        [func = std::forward<F>(f),
         ... captured_args = std::forward<Args>(args)]() mutable -> ReturnType {
          return std::invoke(std::move(func), std::move(captured_args)...);
        });

    std::future<ReturnType> result = task_ptr->get_future();

    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (stop_) {
        throw std::runtime_error("Cannot enqueue on stopped ThreadPool");
      }

      tasks_.emplace([task_ptr]() {
        (*task_ptr)();
      });
    }

    cv_.notify_one();
    return result;
  }

  void shutdown() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (stop_) {
        return;
      }
      stop_ = true;
    }

    cv_.notify_all();

    for (std::thread& worker : workers_) {
      if (worker.joinable()) {
        worker.join();
      }
    }
  }

 private:
  void workerLoop() {
    while (true) {
      std::function<void()> task;

      {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this]() {
          return stop_ || !tasks_.empty();
        });

        if (stop_ && tasks_.empty()) {
          return;
        }

        task = std::move(tasks_.front());
        tasks_.pop();
      }

      task();
    }
  }

  std::vector<std::thread> workers_;
  std::queue<std::function<void()>> tasks_;

  std::mutex mutex_;
  std::condition_variable cv_;
  bool stop_;
};