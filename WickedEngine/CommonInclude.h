#pragma once

// This is a helper include file pasted into all engine headers, try to keep it minimal!
// Do not include engine features in this file!

#include <cassert>
#include <cmath>
#include <cstdint>
#include <type_traits>

#define arraysize(a) (sizeof(a) / sizeof((a)[0]))

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif // NOMINMAX
#endif // _WIN32

#if defined(__GNUC__) || defined(__clang__)
#define __forceinline __attribute__((always_inline)) inline
#define NO_SANITIZE(x) __attribute__((no_sanitize(x)))
#else
#define NO_SANITIZE(x) (void)(x)
#endif // defined(__GNUC__) || defined(__clang__)

#ifdef _MSC_VER
#define WI_DISABLE_DEPRECATED_BEGIN __pragma(warning(push)) __pragma(warning(disable: 4996))
#define WI_DISABLE_DEPRECATED_END __pragma(warning(pop))
#else
#define WI_DISABLE_DEPRECATED_BEGIN _Pragma("GCC diagnostic push") _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#define WI_DISABLE_DEPRECATED_END _Pragma("GCC diagnostic pop")
#endif // _MSC_VER

// Simple common math helpers:

template<typename T>
[[nodiscard]] constexpr T align(T value, T alignment) noexcept
{
	// Fast path for power-of-2 alignments (most common in graphics):
	// Uses bitwise operations: (value + alignment - 1) & ~(alignment - 1)
	// General case for arbitrary alignment:
	// Uses division: ((value + alignment - 1) / alignment) * alignment
	const T mask = alignment - 1;

	if ((alignment & mask) == 0) // power of 2
	{
		return (value + mask) & ~mask;
	}

	// General case - safe for any alignment
	return ((value + alignment - 1) / alignment) * alignment;
}

template<typename T>
[[nodiscard]] constexpr bool is_aligned(T value, T alignment) noexcept
{
	// Fast path for power-of-2 alignments (most common in graphics)
	const T mask = alignment - 1;

	if ((alignment & mask) == 0) // power of 2
	{
		return (value & mask) == 0;
	}

	// General case
	return value % alignment == 0;
}

template <typename T>
[[nodiscard]] constexpr T sqr(T x) noexcept { return x * x; }

template <typename T>
[[nodiscard]] constexpr T pow4(T x) noexcept { return sqr(sqr(x)); }

template <typename T>
[[nodiscard]] constexpr T clamp(T x, T a, T b) noexcept
{
	return x < a ? a : (x > b ? b : x);
}

template <typename T>
[[nodiscard]] constexpr T saturate(T x) noexcept
{
	return clamp(x, T(0), T(1));
}

template <typename T>
[[nodiscard]] constexpr T frac(T x) noexcept
{
	static_assert(std::is_floating_point_v<T>, "frac only supports floating-point types");
	return x - std::floor(x);
}

template <typename T>
[[nodiscard]] constexpr T lerp(T x, T y, T a) noexcept
{
	return x + a * (y - x);
}

template <typename T>
[[nodiscard]] constexpr T inverse_lerp(T value1, T value2, T pos) noexcept
{
	// If value1 == value2, range is zero - return 0 as safe fallback.
	// This handles degenerate cases (zero-size AABBs, flat terrain, etc.)
	return value2 == value1 ? T(0) : ((pos - value1) / (value2 - value1));
}

template <typename T>
[[nodiscard]] constexpr T smoothstep(T edge0, T edge1, T x) noexcept
{
	const T t = saturate((x - edge0) / (edge1 - edge0));
	return t * t * (T(3) - T(2) * t);
}

template <typename Vec4, typename Vec2>
[[nodiscard]] constexpr auto bilinear(Vec4 gather, Vec2 pixel_frac) noexcept
{
	using T = decltype(gather.x);
	const T top_row = lerp(gather.w, gather.z, pixel_frac.x);
	const T bottom_row = lerp(gather.x, gather.y, pixel_frac.x);
	return lerp(top_row, bottom_row, pixel_frac.y);
}

// Stack allocated string utility:
template <unsigned Capacity = 256>
struct StackString
{
	char chars[Capacity] = {};
	unsigned cnt = 0;
	static_assert(Capacity > 1);
	[[nodiscard]] constexpr operator const char* () const noexcept { return chars; }
	[[nodiscard]] constexpr const char* const c_str() const noexcept { return chars; }
	constexpr void push_back(const char* str) noexcept { if (!str) return; while (*str != 0 && (cnt < (Capacity - 1))) { chars[cnt++] = *str; str++; } }
	constexpr void push_back(char c) noexcept { if (cnt < (Capacity - 1)) { chars[cnt++] = c; } }
	[[nodiscard]] constexpr StackString& operator+=(char c) noexcept { push_back(c); return *this; }
	[[nodiscard]] constexpr StackString& operator+=(const char* str) noexcept { push_back(str); return *this; }
	[[nodiscard]] constexpr unsigned size() const noexcept { return cnt; }
	[[nodiscard]] constexpr unsigned length() const noexcept { return cnt; }
	[[nodiscard]] constexpr unsigned capacity() const noexcept { return Capacity; }
	[[nodiscard]] constexpr bool empty() const noexcept { return cnt == 0; }
	constexpr StackString() = default;
	constexpr StackString(const StackString&) = default;
	constexpr StackString(StackString&&) = default;
	constexpr StackString(const char* str) { push_back(str); }
};

// Stack allocated simplified vector container:
template<typename T, unsigned count>
struct StackVector
{
	T items[count] = {};
	unsigned last = 0;
	constexpr void set_size(unsigned size) noexcept { assert(size <= count); last = size; }
	[[nodiscard]] constexpr unsigned size() const noexcept { return last; }
	[[nodiscard]] constexpr bool empty() const noexcept { return last == 0; }
	[[nodiscard]] constexpr unsigned capacity() const noexcept { return count; }
	constexpr void push_back(const T& item) noexcept { assert(last < count); items[last++] = item; }
	constexpr void push_back(T&& item) noexcept { assert(last < count); items[last++] = static_cast<T&&>(item); }
	constexpr void pop_back() noexcept { if (!empty()) items[--last] = {}; }
	constexpr void clear() noexcept { for (unsigned i = 0; i < last; ++i) items[i] = {}; last = 0; }
	[[nodiscard]] constexpr T* data() noexcept { return items; }
	[[nodiscard]] constexpr const T* data() const noexcept { return items; }
	[[nodiscard]] constexpr T* begin() noexcept { return items; }
	[[nodiscard]] constexpr T* end() noexcept { return items + last; }
	constexpr const T& back() const noexcept { assert(!empty()); return items[last - 1]; }
	constexpr T& back() noexcept { assert(!empty()); return items[last - 1]; }
	constexpr const T& front() const noexcept { assert(!empty()); return items[0]; }
	constexpr T& front() noexcept { assert(!empty()); return items[0]; }
	constexpr T& emplace_back() noexcept { push_back({}); return back(); }
	constexpr const T& operator[](unsigned index) const noexcept { assert(index < last); return items[index]; }
	constexpr T& operator[](unsigned index) noexcept { assert(index < last); return items[index]; }
};

// CPU intrinsics:
#if defined(_WIN32)
// Windows, Xbox:
#include <intrin.h>
#if defined(_M_ARM64)
#include <arm64intr.h>
#endif // _M_ARM64
[[nodiscard]] inline long AtomicAnd(volatile long* ptr, long mask) noexcept
{
	return _InterlockedAnd(ptr, mask);
}
[[nodiscard]] inline long long AtomicAnd(volatile long long* ptr, long long mask) noexcept
{
	return _InterlockedAnd64(ptr, mask);
}
[[nodiscard]] inline long AtomicOr(volatile long* ptr, long mask) noexcept
{
	return _InterlockedOr(ptr, mask);
}
[[nodiscard]] inline long long AtomicOr(volatile long long* ptr, long long mask) noexcept
{
	return _InterlockedOr64(ptr, mask);
}
[[nodiscard]] inline long AtomicXor(volatile long* ptr, long mask) noexcept
{
	return _InterlockedXor(ptr, mask);
}
[[nodiscard]] inline long long AtomicXor(volatile long long* ptr, long long mask) noexcept
{
	return _InterlockedXor64(ptr, mask);
}
[[nodiscard]] inline long AtomicAdd(volatile long* ptr, long val) noexcept
{
	return _InterlockedExchangeAdd(ptr, val);
}
[[nodiscard]] inline long long AtomicAdd(volatile long long* ptr, long long val) noexcept
{
	return _InterlockedExchangeAdd64(ptr, val);
}
[[nodiscard]] inline unsigned int countbits(unsigned int value) noexcept
{
	return __popcnt(value);
}
[[nodiscard]] inline unsigned long countbits(unsigned long value) noexcept
{
	return (unsigned long)__popcnt64((unsigned long long)value);
}
[[nodiscard]] inline unsigned long long countbits(unsigned long long value) noexcept
{
	return __popcnt64(value);
}
[[nodiscard]] inline unsigned int firstbithigh(unsigned int value) noexcept
{
	unsigned long bit_index;
	if (_BitScanReverse(&bit_index, (unsigned long)value))
	{
		return 31u - (unsigned int)bit_index;
	}
	return 32u;
}
[[nodiscard]] inline unsigned long firstbithigh(unsigned long value) noexcept
{
	unsigned long bit_index;
	if (_BitScanReverse64(&bit_index, (unsigned long long)value))
	{
		return 31ul - bit_index;
	}
	return 32ul;
}
[[nodiscard]] inline unsigned long long firstbithigh(unsigned long long value) noexcept
{
	unsigned long bit_index;
	if (_BitScanReverse64(&bit_index, value))
	{
		return 63ull - bit_index;
	}
	return 64ull;
}
[[nodiscard]] inline unsigned int firstbitlow(unsigned int value) noexcept
{
	unsigned long bit_index;
	if (_BitScanForward(&bit_index, value))
	{
		return (unsigned int)bit_index;
	}
	return 32u;
}
[[nodiscard]] inline unsigned long firstbitlow(unsigned long value) noexcept
{
	unsigned long bit_index;
	if (_BitScanForward64(&bit_index, (unsigned long long)value))
	{
		return bit_index;
	}
	return 32ul;
}
[[nodiscard]] inline unsigned long firstbitlow(unsigned long long value) noexcept
{
	unsigned long bit_index;
	if (_BitScanForward64(&bit_index, value))
	{
		return bit_index;
	}
	return 64ul;
}
#else
// Linux, PlayStation:
[[nodiscard]] inline long AtomicAnd(volatile long* ptr, long mask) noexcept
{
	return __atomic_fetch_and(ptr, mask, __ATOMIC_SEQ_CST);
}
[[nodiscard]] inline long long AtomicAnd(volatile long long* ptr, long long mask) noexcept
{
	return __atomic_fetch_and(ptr, mask, __ATOMIC_SEQ_CST);
}
[[nodiscard]] inline long AtomicOr(volatile long* ptr, long mask) noexcept
{
	return __atomic_fetch_or(ptr, mask, __ATOMIC_SEQ_CST);
}
[[nodiscard]] inline long long AtomicOr(volatile long long* ptr, long long mask) noexcept
{
	return __atomic_fetch_or(ptr, mask, __ATOMIC_SEQ_CST);
}
[[nodiscard]] inline long AtomicXor(volatile long* ptr, long mask) noexcept
{
	return __atomic_fetch_xor(ptr, mask, __ATOMIC_SEQ_CST);
}
[[nodiscard]] inline long long AtomicXor(volatile long long* ptr, long long mask) noexcept
{
	return __atomic_fetch_xor(ptr, mask, __ATOMIC_SEQ_CST);
}
[[nodiscard]] inline long AtomicAdd(volatile long* ptr, long val) noexcept
{
	return __atomic_fetch_add(ptr, val, __ATOMIC_SEQ_CST);
}
[[nodiscard]] inline long long AtomicAdd(volatile long long* ptr, long long val) noexcept
{
	return __atomic_fetch_add(ptr, val, __ATOMIC_SEQ_CST);
}
[[nodiscard]] inline unsigned int countbits(unsigned int value) noexcept
{
	return __builtin_popcount(value);
}
[[nodiscard]] inline unsigned long countbits(unsigned long value) noexcept
{
	return __builtin_popcountl(value);
}
[[nodiscard]] inline unsigned long long countbits(unsigned long long value) noexcept
{
	return __builtin_popcountll(value);
}
[[nodiscard]] inline unsigned int firstbithigh(unsigned int value) noexcept
{
	if (value == 0)
	{
		return 32u;
	}
	return __builtin_clz(value);
}
[[nodiscard]] inline unsigned long firstbithigh(unsigned long value) noexcept
{
	if (value == 0)
	{
		return sizeof(unsigned long) * 8;
	}
	return __builtin_clzl(value);
}
[[nodiscard]] inline unsigned long long firstbithigh(unsigned long long value) noexcept
{
	if (value == 0)
	{
		return 64ull;
	}
	return __builtin_clzll(value);
}
[[nodiscard]] inline unsigned int firstbitlow(unsigned int value) noexcept
{
	if (value == 0)
	{
		return 32u;
	}
	return __builtin_ctz(value);
}
[[nodiscard]] inline unsigned long firstbitlow(unsigned long value) noexcept
{
	if (value == 0)
	{
		return sizeof(unsigned long) * 8;
	}
	return __builtin_ctzl(value);
}
[[nodiscard]] inline unsigned long long firstbitlow(unsigned long long value) noexcept
{
	if (value == 0)
	{
		return 64ull;
	}
	return __builtin_ctzll(value);
}
#endif // _WIN32

[[nodiscard]] inline long AtomicLoad(const volatile long* ptr) noexcept
{
	return AtomicOr((volatile long*)ptr, 0);
}
[[nodiscard]] inline long long AtomicLoad(const volatile long long* ptr) noexcept
{
	return AtomicOr((volatile long long*)ptr, 0);
}

// Enable enum flags:
//	https://www.justsoftwaresolutions.co.uk/cplusplus/using-enum-classes-as-bitfields.html
template<typename E>
struct enable_bitmask_operators : std::false_type {};

template<typename E>
using EnableIfBitmaskOps = std::enable_if_t<enable_bitmask_operators<E>::value, int>;

template<typename E, EnableIfBitmaskOps<E> = 0>
constexpr E operator|(E lhs, E rhs) noexcept
{
	using U = std::underlying_type_t<E>;
	return static_cast<E>(static_cast<U>(lhs) | static_cast<U>(rhs));
}
template<typename E, EnableIfBitmaskOps<E> = 0>
constexpr E operator|=(E& lhs, E rhs) noexcept
{
	using U = std::underlying_type_t<E>;
	lhs = static_cast<E>(static_cast<U>(lhs) | static_cast<U>(rhs));
	return lhs;
}
template<typename E, EnableIfBitmaskOps<E> = 0>
constexpr E operator&(E lhs, E rhs) noexcept
{
	using U = std::underlying_type_t<E>;
	return static_cast<E>(static_cast<U>(lhs) & static_cast<U>(rhs));
}
template<typename E, EnableIfBitmaskOps<E> = 0>
constexpr E operator&=(E& lhs, E rhs) noexcept
{
	using U = std::underlying_type_t<E>;
	lhs = static_cast<E>(static_cast<U>(lhs) & static_cast<U>(rhs));
	return lhs;
}
template<typename E, EnableIfBitmaskOps<E> = 0>
constexpr E operator~(E rhs) noexcept
{
	using U = std::underlying_type_t<E>;
	rhs = static_cast<E>(~static_cast<U>(rhs));
	return rhs;
}
template<typename E, EnableIfBitmaskOps<E> = 0>
[[nodiscard]] constexpr bool has_flag(E lhs, E rhs) noexcept
{
	return (lhs & rhs) == rhs;
}
template<typename T, typename U>
constexpr void set_flag(T& flags, U flag, bool set) noexcept
{
	if (set)
	{
		flags |= static_cast<T>(flag);
	}
	else
	{
		flags &= ~static_cast<T>(flag);
	}
}
// Extract file name from a path at compile-time
[[nodiscard]] constexpr const char* relative_path(const char* path) noexcept
{
	const char* startPosition = path;
	for (const char* currentCharacter = path; *currentCharacter != '\0'; ++currentCharacter)
	{
		if (*currentCharacter == '\\' || *currentCharacter == '/')
		{
			startPosition = currentCharacter;
		}
	}

	if (startPosition != path)
	{
		++startPosition;
	}

	return startPosition;
}

[[nodiscard]] constexpr auto relative_path_storage(const char* path) noexcept
{
	StackString ret;
	ret.push_back(relative_path(path));
	return ret;
}

// Convert string literal to StackString at compile-time
[[nodiscard]] constexpr auto to_stack_string(const char* str) noexcept
{
	StackString ret;
	ret.push_back(str);
	return ret;
}
