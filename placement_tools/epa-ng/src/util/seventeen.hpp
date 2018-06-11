#pragma once

#include <type_traits>

namespace seventeen {

template <typename T, typename F>
auto static_if(std::true_type, T t, F) { return t; }

template <typename T, typename F>
auto static_if(std::false_type, T, F f) { return f; }

template <bool B, typename T, typename F>
auto static_if(T t, F f) { return static_if(std::integral_constant<bool, B>{}, t, f); }

template <bool B, typename T>
auto static_if(T t) { return static_if(std::integral_constant<bool, B>{}, t, [](auto&&...){}); }

}