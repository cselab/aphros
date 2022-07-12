#pragma once

#include <memory>
#include <utility>

template<typename T, typename... Args>
std::unique_ptr<T> MakeUnique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
