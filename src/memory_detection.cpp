#include <cstdlib>
#include <chrono>
#include <thread>
#include <iostream>
#include <fstream>
#include <csignal>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>

// 每次采样间隔（毫秒）
const int SAMPLE_INTERVAL_MS = 50;

int main() {
    // Step 1: Fork 子进程执行 ./GPU_projection
    pid_t pid = fork();

    if (pid == 0) {
        // 子进程：运行你的 CUDA 程序
        execl("./GPU_projection", "./GPU_projection", nullptr);
        // 如果运行失败
        std::cerr << "Failed to launch ./GPU_projection" << std::endl;
        std::exit(1);
    }

    // Step 2: 父进程：开始采样 nvidia-smi
    std::cout << "Started ./GPU_projection with PID " << pid << std::endl;
    std::ofstream log("mem_log.txt");

    if (!log) {
        std::cerr << "Failed to open mem_log.txt for writing" << std::endl;
        return 1;
    }

    while (true) {
        // 检查 ./GPU_projection 是否还在运行
        int status;
        pid_t result = waitpid(pid, &status, WNOHANG);

        if (result == 0) {
            // 子进程还活着，继续采样
            FILE* pipe = popen("nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader,nounits", "r");
            if (pipe) {
                char buffer[256];
                while (fgets(buffer, sizeof(buffer), pipe)) {
                    // 仅记录目标子进程的 GPU 显存
                    std::string line(buffer);
                    if (line.find(std::to_string(pid)) != std::string::npos) {
                        log << line;
                    }
                }
                pclose(pipe);
            }

            std::this_thread::sleep_for(std::chrono::milliseconds(SAMPLE_INTERVAL_MS));
        } else {
            // 子进程结束
            std::cout << "\n[Info] ./GPU_projection has exited. Stopping memory logging.\n";
            break;
        }
    }

    log.close();
    return 0;
}
