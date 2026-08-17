/**
 * WickedEngine Common Include Header
 *
 * This header provides core utilities, math helpers, containers, atomic
 * operations, and bit manipulation functions used throughout the WickedEngine.
 * It is designed to be minimal and included in every engine header.
 *
 * WARNING: Do not include engine features in this file!
 */

#pragma once

#include <cassert>
#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <numeric>
#include <type_traits>

/**
 * Array size helper macro.
 *
 * Returns the number of elements in a statically-sized array.
 *
 * @param a - Array expression (not a pointer).
 *
 * @return Number of elements in the array.
 */
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
#define WI_DISABLE_DEPRECATED_BEGIN \
	__pragma(warning(push)) \
	__pragma(warning(disable: 4996))
#define WI_DISABLE_DEPRECATED_END \
	__pragma(warning(pop))
#else
#define WI_DISABLE_DEPRECATED_BEGIN \
	_Pragma("GCC diagnostic push") \
	_Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#define WI_DISABLE_DEPRECATED_END \
	_Pragma("GCC diagnostic pop")
#endif // _MSC_VER

/*
################################################################################
Math Helpers
################################################################################
*/

/**
 * Aligns a value up to the specified alignment.
 *
 * Supports both power-of-2 and arbitrary alignments.
 * For power-of-2 alignments, uses fast bitwise operations.
 * For arbitrary alignments, uses division-based approach.
 *
 * \[
 * \text{align}(v, a) = \lceil \frac{v}{a} \rceil \times a
 * \]
 *
 * @param[in] value - Value to align up.
 * @param[in] alignment - Alignment boundary (must be > 0).
 *
 * @return The smallest multiple of `alignment` that is >= `value`.
 *
 * @note Fast path for power-of-2 alignments uses bitwise operations.
 * @pre `alignment != 0`
 */
template<typename T>
[[nodiscard]] constexpr T align(T value, T alignment) noexcept
{
	static_assert(std::is_unsigned_v<T>, "align requires unsigned integer type");
	assert(alignment != 0);

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

/**
 * Checks if a value is aligned to the specified alignment.
 *
 * Supports both power-of-2 and arbitrary alignments.
 *
 * @param[in] value - Value to check.
 * @param[in] alignment - Alignment boundary (must be > 0).
 *
 * @return true if `value` is a multiple of `alignment`, false otherwise.
 *
 * @pre `alignment != 0`
 */
template<typename T>
[[nodiscard]] constexpr bool is_aligned(T value, T alignment) noexcept
{
	static_assert(std::is_unsigned_v<T>, "is_aligned requires unsigned integer type");
	assert(alignment != 0);

	// Fast path for power-of-2 alignments (most common in graphics)
	const T mask = alignment - 1;

	if ((alignment & mask) == 0) // power of 2
	{
		return (value & mask) == 0;
	}

	// General case
	return value % alignment == 0;
}

/**
 * Computes the square of a value.
 *
 * @param[in] x - Input value.
 *
 * @return \f$ x^2 \f$
 */
template <typename T>
[[nodiscard]] constexpr T sqr(T x) noexcept
{
	return x * x;
}

/**
 * Computes the fourth power of a value.
 *
 * @param[in] x - Input value.
 *
 * @return \f$ x^4 \f$
 */
template <typename T>
[[nodiscard]] constexpr T pow4(T x) noexcept
{
	return sqr(sqr(x));
}

/**
 * Clamps a value between two bounds.
 *
 * @param[in] x - Value to clamp.
 * @param[in] a - Lower bound.
 * @param[in] b - Upper bound.
 *
 * @return Value clamped to range [a, b].
 */
template <typename T>
[[nodiscard]] constexpr T clamp(T x, T a, T b) noexcept
{
	return std::clamp(x, a, b);
}

/**
 * Clamps a value to the [0, 1] range.
 *
 * @param[in] x - Value to saturate.
 *
 * @return Value clamped to [0, 1].
 */
template <typename T>
[[nodiscard]] constexpr T saturate(T x) noexcept
{
	return clamp(x, T(0), T(1));
}

/**
 * Returns the fractional part of a floating-point value.
 *
 * Returns value in range [0, 1).
 *
 * @param[in] x - Input value (must be floating-point).
 *
 * @return Fractional part of x.
 */
template <typename T>
[[nodiscard]] constexpr T frac(T x) noexcept
{
	static_assert(
		std::is_floating_point_v<T>,
		"frac only supports floating-point types"
	);
	T intpart;

	T f = std::modf(x, &intpart);

	// std::modf returns fractional part with same sign as x.
	// For negative x, we need to add 1 to get [0, 1) range.
	return f < 0 ? f + T(1) : f;
}

/**
 * Linear interpolation between two values.
 *
 * Uses numerically stable formula: \f$ x + a \times (y - x) \f$
 *
 * @param[in] x - Start value.
 * @param[in] y - End value.
 * @param[in] a - Interpolation factor in [0, 1].
 *
 * @return Interpolated value.
 */
template <typename T>
[[nodiscard]] constexpr T lerp(T x, T y, T a) noexcept
{
	return std::lerp(x, y, a);
}

/**
 * Computes the interpolation factor between two values.
 *
 * Returns 0 if value1 == value2 (zero range).
 *
 * @param[in] value1 - Start value.
 * @param[in] value2 - End value.
 * @param[in] pos - Position to find factor for.
 *
 * @return Factor a such that lerp(value1, value2, a) == pos.
 */
template <typename T>
[[nodiscard]] constexpr T inverse_lerp(T value1, T value2, T pos) noexcept
{
	// If value1 == value2, range is zero - return 0 as safe fallback.
	// This handles degenerate cases (zero-size AABBs, flat terrain, etc.)
	return value2 == value1 ? T(0) : ((pos - value1) / (value2 - value1));
}

/**
 * Smooth Hermite interpolation between 0 and 1.
 *
 * \[
 * f(t) = t^2(3 - 2t)
 * \]
 *
 * @param[in] edge0 - Lower edge.
 * @param[in] edge1 - Upper edge.
 * @param[in] x - Input value.
 *
 * @return Smooth interpolated value in [0, 1].
 */
template <typename T>
[[nodiscard]] constexpr T smoothstep(T edge0, T edge1, T x) noexcept
{
	const T t = saturate((x - edge0) / (edge1 - edge0));

	return t * t * (T(3) - T(2) * t);
}
 
/**
 * Type trait to detect if a type has .x, .y arithmetic members (Vec2-like).
 */
template <typename T, typename = void>
struct is_vec2_like : std::false_type {};
 
template <typename T>
struct is_vec2_like<T,
    std::void_t<
        decltype(std::declval<T>().x),
        decltype(std::declval<T>().y)
    >
> : std::integral_constant<bool,
    std::is_arithmetic_v<decltype(std::declval<T>().x)> &&
    std::is_arithmetic_v<decltype(std::declval<T>().y)>
> {};
 
/**
 * Type trait to detect if a type has .x, .y, .z, .w arithmetic members
 * (Vec4-like).
 */
template <typename T, typename = void>
struct is_vec4_like : std::false_type {};
 
template <typename T>
struct is_vec4_like<T,
    std::void_t<
        decltype(std::declval<T>().x),
        decltype(std::declval<T>().y),
        decltype(std::declval<T>().z),
        decltype(std::declval<T>().w)
    >
> : std::integral_constant<bool,
    std::is_arithmetic_v<decltype(std::declval<T>().x)> &&
    std::is_arithmetic_v<decltype(std::declval<T>().y)> &&
    std::is_arithmetic_v<decltype(std::declval<T>().z)> &&
    std::is_arithmetic_v<decltype(std::declval<T>().w)>
> {};
 
/**
 * Bilinear interpolation on a 2x2 grid.
 *
 * Interpolates between four corner values using a 2D fractional position.
 *
 * @param[in] gather - 4-component vector containing corner values
 *                     [bottom-left, bottom-right, top-left, top-right].
 * @param[in] pixel_frac - 2D fractional position in [0, 1] x [0, 1].
 *
 * @return Bilinearly interpolated value.
 *
 * @note Vec4 must have .x, .y, .z, .w arithmetic members.
 * @note Vec2 must have .x, .y arithmetic members.
 */
template <typename Vec4, typename Vec2,
    std::enable_if_t<is_vec4_like<Vec4>::value && is_vec2_like<Vec2>::value, int> = 0>
[[nodiscard]] constexpr auto bilinear(Vec4 gather, Vec2 pixel_frac) noexcept
    -> decltype(gather.x)
{
	using T = decltype(gather.x);

	const T top_row    = lerp(gather.w, gather.z, pixel_frac.x);
	const T bottom_row = lerp(gather.x, gather.y, pixel_frac.x);

	return lerp(top_row, bottom_row, pixel_frac.y);
}

/**
 * Stack-allocated string utility with fixed capacity.
 *
 * A lightweight string container that stores characters on the stack.
 * Capacity is fixed at compile-time via template parameter.
 * Not dynamically resizable - use with care for large strings.
 *
 * @tparam Capacity - Maximum number of characters (including null terminator).
 *         Must be > 1.
 *
 * Example usage:
 * ```cpp
 * StackString<64> str = "hello";
 * str += " world";
 * str += '!';
 * printf("%s\n", str.c_str()); // "hello world!"
 * ```
 */
template <std::size_t Capacity = 256>
struct StackString
{
	/** Character storage buffer (including null terminator). */
	char chars[Capacity] = {};

	/**
	 * Current number of characters in the string (excludes null terminator).
	 */
	std::size_t cnt = 0;

	static_assert(Capacity > 1);
	
	/**
	 * Implicit conversion to C-string.
	 *
	 * @return Null-terminated C-string pointer.
	 */
	[[nodiscard]] constexpr operator const char* () const noexcept
	{
		return chars;
	}
	
	/**
	 * Returns the C-string pointer.
	 *
	 * @return Null-terminated C-string pointer.
	 */
	[[nodiscard]] constexpr const char* const c_str() const noexcept
	{
		return chars;
	}
	
	/**
	 * Appends a null-terminated string.
	 *
	 * @param[in] str - Null-terminated string to append. nullptr is ignored.
	 */
	constexpr void push_back(const char* str) noexcept
	{
		if (!str)
		{
			return;
		}

		while (*str != 0 && (cnt < (Capacity - 1)))
		{
			chars[cnt++] = *str;
			str++;
		}
	}
	
	/**
	 * Appends a single character.
	 *
	 * @param[in] c - Character to append.
	 */
	constexpr void push_back(char c) noexcept
	{
		if (cnt < (Capacity - 1))
		{
			chars[cnt++] = c;
		}
	}
	
	/**
	 * Appends a single character (operator overload).
	 *
	 * @param[in] c - Character to append.
	 *
	 * @return Reference to this StackString for chaining.
	 */
	[[nodiscard]] constexpr StackString& operator+=(char c) noexcept
	{
		push_back(c);

		return *this;
	}
	
	/**
	 * Appends a null-terminated string (operator overload).
	 *
	 * @param[in] str - Null-terminated string to append.
	 *
	 * @return Reference to this StackString for chaining.
	 */
	[[nodiscard]] constexpr StackString& operator+=(const char* str) noexcept
	{
		push_back(str);

		return *this;
	}
	
	/**
	 * Returns the current string length (number of characters).
	 *
	 * @return Number of characters in the string.
	 */
[[nodiscard]] constexpr std::size_t size() const noexcept
{
	return cnt;
}
	
/**
 * Returns the current string length (alias for size()).
 *
 * @return Number of characters in the string.
 */
[[nodiscard]] constexpr std::size_t length() const noexcept
{
	return cnt;
}
	
/**
 * Returns the compile-time capacity of the string.
 *
 * @return Maximum number of characters (including null terminator).
 */
[[nodiscard]] constexpr std::size_t capacity() const noexcept
{
	return Capacity;
}
	
	/**
	 * Checks if the string is empty.
	 *
	 * @return true if string length is 0, false otherwise.
	 */
	[[nodiscard]] constexpr bool empty() const noexcept
	{
		return cnt == 0;
	}
	
	/**
	 * Default constructor.
	 *
	 * Creates an empty string.
	 */
	constexpr StackString() noexcept = default;

	/**
	 * Copy constructor.
	 *
	 * Performs memberwise copy of the string data.
	 */
	constexpr StackString(const StackString&) noexcept = default;

	/**
	 * Move constructor.
	 *
	 * Performs memberwise move of the string data.
	 */
	constexpr StackString(StackString&&) noexcept = default;
	
	/**
	 * Constructs from a C-string.
	 *
	 * @param[in] str - Null-terminated string to copy.
	 */
	constexpr StackString(const char* str) noexcept
	{
		push_back(str);
	}
};

/**
 * Stack-allocated simplified vector container with fixed capacity.
 *
 * A lightweight vector container that stores elements on the stack.
 * Capacity is fixed at compile-time via template parameter.
 * Provides basic vector operations with bounds checking via asserts.
 *
 * @tparam T - Element type.
 * @tparam count - Maximum number of elements. Must be > 0.
 *
 * Example usage:
 * ```cpp
 * StackVector<int, 32> vec;
 * vec.push_back(1);
 * vec.push_back(2);
 * vec.push_back(3);
 * for (int v : vec) { printf("%d ", v); } // 1 2 3
 * ```
 */
template<typename T, std::size_t count>
struct StackVector
{
	/** Element storage array (fixed compile-time capacity). */
	T items[count] = {};

	/** Current number of elements in the vector (logical size). */
	std::size_t last = 0;
	
	/**
	 * Sets the logical size of the vector.
	 *
	 * @param[in] size - New logical size. Must not exceed compile-time
	 *                   capacity.
	 */
	constexpr void set_size(std::size_t size) noexcept
	{
		assert(size <= count);

		last = size;
	}
	
	/**
	 * Returns the current number of elements.
	 *
	 * @return Number of elements in the vector.
	 */
	[[nodiscard]] constexpr std::size_t size() const noexcept
	{
		return last;
	}
	
	/**
	 * Checks if the vector is empty.
	 *
	 * @return true if size is 0, false otherwise.
	 */
	[[nodiscard]] constexpr bool empty() const noexcept
	{
		return last == 0;
	}
	
	/**
	 * Returns the compile-time capacity of the vector.
	 *
	 * @return Maximum number of elements.
	 */
	[[nodiscard]] constexpr std::size_t capacity() const noexcept
	{
		return count;
	}
	
	/**
	 * Appends an element by copy.
	 *
	 * @param[in] item - Element to append.
	 *
	 * @note Asserts if vector is full.
	 */
	constexpr void push_back(const T& item) noexcept
	{
		assert(last < count);

		items[last++] = item;
	}
	
	/**
	 * Appends an element by move.
	 *
	 * @param[in] item - Element to append (rvalue reference).
	 *
	 * @note Asserts if vector is full.
	 */
	constexpr void push_back(T&& item) noexcept
	{
		assert(last < count);

		items[last++] = static_cast<T&&>(item);
	}
	
	/**
	 * Removes the last element.
	 *
	 * @note No-op if vector is empty.
	 */
	constexpr void pop_back() noexcept
	{
		if (!empty())
		{
			items[--last] = {};
		}
	}
	
	/**
	 * Clears all elements (sets size to 0).
	 */
	constexpr void clear() noexcept
	{
		for (std::size_t i = 0; i < last; ++i)
		{
			items[i] = {};
		}

		last = 0;
	}
	
	/**
	 * Returns pointer to the first element.
	 *
	 * @return Pointer to the first element.
	 */
	[[nodiscard]] constexpr T* data() noexcept
	{
		return items;
	}
	
	/**
	 * Returns const pointer to the first element.
	 *
	 * @return Const pointer to the first element.
	 */
	[[nodiscard]] constexpr const T* data() const noexcept
	{
		return items;
	}
	
	/**
	 * Returns iterator to the first element.
	 *
	 * @return Pointer to the first element.
	 */
	[[nodiscard]] constexpr T* begin() noexcept
	{
		return items;
	}
	
	/**
	 * Returns iterator past the last element.
	 *
	 * @return Pointer past the last element.
	 */
	[[nodiscard]] constexpr T* end() noexcept
	{
		return items + last;
	}
	
	/**
	 * Returns reference to the last element.
	 *
	 * @return Reference to the last element.
	 *
	 * @note Asserts if vector is empty.
	 */
	constexpr const T& back() const noexcept
	{
		assert(!empty());

		return items[last - 1];
	}
	
	/**
	 * Returns mutable reference to the last element.
	 *
	 * @return Reference to the last element.
	 *
	 * @note Asserts if vector is empty.
	 */
	constexpr T& back() noexcept
	{
		assert(!empty());

		return items[last - 1];
	}
	
	/**
	 * Returns reference to the first element.
	 *
	 * @return Reference to the first element.
	 *
	 * @note Asserts if vector is empty.
	 */
	constexpr const T& front() const noexcept
	{
		assert(!empty());

		return items[0];
	}
	
	/**
	 * Returns mutable reference to the first element.
	 *
	 * @return Reference to the first element.
	 *
	 * @note Asserts if vector is empty.
	 */
	constexpr T& front() noexcept
	{
		assert(!empty());

		return items[0];
	}
	
	/**
	 * Constructs an element in-place at the end.
	 *
	 * @return Reference to the newly constructed element.
	 */
	constexpr T& emplace_back() noexcept
	{
		push_back({});

		return back();
	}
	
	/**
	 * Returns const reference to element at index.
	 *
	 * @param[in] index - Element index.
	 *
	 * @return Const reference to element.
	 *
	 * @note Asserts if index >= size().
	 */
	constexpr const T& operator[](std::size_t index) const noexcept
	{
		assert(index < last);

		return items[index];
	}
	
	/**
	 * Returns mutable reference to element at index.
	 *
	 * @param[in] index - Element index.
	 *
	 * @return Reference to element.
	 *
	 * @note Asserts if index >= size().
	 */
	constexpr T& operator[](std::size_t index) noexcept
	{
		assert(index < last);

		return items[index];
	}
	
	/**
	 * Returns const reference to element at index (bounds-checked alias).
	 *
	 * @param[in] index - Element index.
	 *
	 * @return Const reference to element.
	 *
	 * @note Asserts if index >= size().
	 */
	constexpr const T& at(std::size_t index) const noexcept
	{
		assert(index < last);

		return items[index];
	}
	
	/**
	 * Returns mutable reference to element at index (bounds-checked alias).
	 *
	 * @param[in] index - Element index.
	 *
	 * @return Reference to element.
	 *
	 * @note Asserts if index >= size().
	 */
	constexpr T& at(std::size_t index) noexcept
	{
		assert(index < last);

		return items[index];
	}
	
	};

/*
################################################################################
Atomic Operations
################################################################################
*/

/**
 * Atomically performs bitwise AND on an integer.
 *
 * @param[in,out] ptr - Pointer to the value.
 * @param[in] mask - Bitmask to AND with.
 *
 * @return The original value before the operation.
 */
template <typename T, typename U>
[[nodiscard]] inline T AtomicAnd(volatile T* ptr, U mask) noexcept
{
	static_assert(std::is_integral_v<T> && (sizeof(T) == 4 || sizeof(T) == 8),
		"AtomicAnd requires 32-bit or 64-bit integral type");
	static_assert(std::is_integral_v<U>, "Mask must be integral type");

	std::atomic_ref<T> ref(*const_cast<T*>(ptr));

	return ref.fetch_and(static_cast<T>(mask), std::memory_order_acq_rel);
}

/**
 * Atomically performs bitwise OR on an integer.
 *
 * @param[in,out] ptr - Pointer to the value.
 * @param[in] mask - Bitmask to OR with.
 *
 * @return The original value before the operation.
 */
template <typename T, typename U>
[[nodiscard]] inline T AtomicOr(volatile T* ptr, U mask) noexcept
{
	static_assert(std::is_integral_v<T> && (sizeof(T) == 4 || sizeof(T) == 8),
		"AtomicOr requires 32-bit or 64-bit integral type");
	static_assert(std::is_integral_v<U>, "Mask must be integral type");

	std::atomic_ref<T> ref(*const_cast<T*>(ptr));

	return ref.fetch_or(static_cast<T>(mask), std::memory_order_acq_rel);
}

/**
 * Atomically performs bitwise XOR on an integer.
 *
 * @param[in,out] ptr - Pointer to the value.
 * @param[in] mask - Bitmask to XOR with.
 *
 * @return The original value before the operation.
 */
template <typename T, typename U>
[[nodiscard]] inline T AtomicXor(volatile T* ptr, U mask) noexcept
{
	static_assert(std::is_integral_v<T> && (sizeof(T) == 4 || sizeof(T) == 8),
		"AtomicXor requires 32-bit or 64-bit integral type");
	static_assert(std::is_integral_v<U>, "Mask must be integral type");

	std::atomic_ref<T> ref(*const_cast<T*>(ptr));

	return ref.fetch_xor(static_cast<T>(mask), std::memory_order_acq_rel);
}

/**
 * Atomically adds a value to an integer.
 *
 * @param[in,out] ptr - Pointer to the value.
 * @param[in] val - Value to add.
 *
 * @return The original value before the addition.
 */
template <typename T, typename U>
[[nodiscard]] inline T AtomicAdd(volatile T* ptr, U val) noexcept
{
	static_assert(std::is_integral_v<T> && (sizeof(T) == 4 || sizeof(T) == 8),
		"AtomicAdd requires 32-bit or 64-bit integral type");
	static_assert(std::is_integral_v<U>, "Value must be integral type");

	std::atomic_ref<T> ref(*const_cast<T*>(ptr));

	return ref.fetch_add(static_cast<T>(val), std::memory_order_acq_rel);
}

/**
 * Atomically loads an integer.
 *
 * @param[in] ptr - Pointer to the value.
 *
 * @return The loaded value.
 */
template <typename T>
[[nodiscard]] inline T AtomicLoad(const volatile T* ptr) noexcept
{
	static_assert(std::is_integral_v<T> && (sizeof(T) == 4 || sizeof(T) == 8),
		"AtomicLoad requires 32-bit or 64-bit integral type");

	// Cast away volatile for std::atomic_ref (which doesn't accept volatile)
	std::atomic_ref<const T> ref(*const_cast<const T*>(ptr));

	return ref.load(std::memory_order_acquire);
}

/*
################################################################################
Bit Manipulation (Linux/PlayStation)
################################################################################
*/

/**
 * Counts the number of set bits (population count) in an unsigned integer.
 *
 * @param[in] value - Unsigned integer value.
 *
 * @return Number of set bits (0 to bit width of type).
 */
template <typename T>
[[nodiscard]] inline T countbits(T value) noexcept
{
	static_assert(std::is_unsigned_v<T>, "countbits requires unsigned integer type");

	return static_cast<T>(std::popcount(value));
}

/**
 * Finds the index of the most significant set bit (0-based from MSB).
 *
 * Returns bit width if value is 0.
 *
 * @param[in] value - Unsigned integer value.
 *
 * @return Index of most significant set bit, or bit width if value is 0.
 *
 * @note Returns bit width for 0 input to distinguish from bit 0 being set.
 */
template <typename T>
[[nodiscard]] inline T firstbithigh(T value) noexcept
{
	static_assert(std::is_unsigned_v<T>, "firstbithigh requires unsigned integer type");

	if (value == 0)
	{
		return static_cast<T>(sizeof(T) * 8);
	}

	return static_cast<T>((sizeof(T) * 8 - 1) - std::countl_zero(value));
}

/**
 * Finds the index of the least significant set bit (0-based from LSB).
 *
 * Returns bit width if value is 0.
 *
 * @param[in] value - Unsigned integer value.
 *
 * @return Index of least significant set bit, or bit width if value is 0.
 */
template <typename T>
[[nodiscard]] inline T firstbitlow(T value) noexcept
{
	static_assert(std::is_unsigned_v<T>, "firstbitlow requires unsigned integer type");

	if (value == 0)
	{
		return static_cast<T>(sizeof(T) * 8);
	}

	return static_cast<T>(std::countr_zero(value));
}

/*
################################################################################
Bitmask Operators
################################################################################
*/

// Enable enum flags:
//	https://www.justsoftwaresolutions.co.uk/cplusplus/using-enum-classes-as-bitfields.html

/**
 * Trait to enable bitmask operators for an enum type.
 *
 * Specialize this to `std::true_type` for enums that should support bitmask
 * operations.
 *
 * @tparam E - Enum type.
 */
template<typename E>
struct enable_bitmask_operators : std::false_type {};

/**
 * SFINAE helper for enabling bitmask operators.
 *
 * @tparam E - Enum type.
 */
template<typename E>
using EnableIfBitmaskOps = std::enable_if_t<
	enable_bitmask_operators<E>::value, int
>;

/**
 * Bitwise OR for bitmask enums.
 *
 * @param[in] lhs - Left-hand side.
 * @param[in] rhs - Right-hand side.
 *
 * @return Result of bitwise OR.
 */
template<typename E, EnableIfBitmaskOps<E> = 0>
constexpr E operator|(E lhs, E rhs) noexcept
{
	using U = std::underlying_type_t<E>;

	return static_cast<E>(static_cast<U>(lhs) | static_cast<U>(rhs));
}

/**
 * Bitwise OR-assign for bitmask enums.
 *
 * @param[in,out] lhs - Left-hand side (modified).
 * @param[in] rhs - Right-hand side.
 *
 * @return Reference to modified lhs.
 */
template<typename E, EnableIfBitmaskOps<E> = 0>
constexpr E operator|=(E& lhs, E rhs) noexcept
{
	using U = std::underlying_type_t<E>;

	lhs = static_cast<E>(static_cast<U>(lhs) | static_cast<U>(rhs));

	return lhs;
}

/**
 * Bitwise AND for bitmask enums.
 *
 * @param[in] lhs - Left-hand side.
 * @param[in] rhs - Right-hand side.
 *
 * @return Result of bitwise AND.
 */
template<typename E, EnableIfBitmaskOps<E> = 0>
constexpr E operator&(E lhs, E rhs) noexcept
{
	using U = std::underlying_type_t<E>;

	return static_cast<E>(static_cast<U>(lhs) & static_cast<U>(rhs));
}

/**
 * Bitwise AND-assign for bitmask enums.
 *
 * @param[in,out] lhs - Left-hand side (modified).
 * @param[in] rhs - Right-hand side.
 *
 * @return Reference to modified lhs.
 */
template<typename E, EnableIfBitmaskOps<E> = 0>
constexpr E operator&=(E& lhs, E rhs) noexcept
{
	using U = std::underlying_type_t<E>;

	lhs = static_cast<E>(static_cast<U>(lhs) & static_cast<U>(rhs));

	return lhs;
}

/**
 * Bitwise NOT for bitmask enums.
 *
 * @param[in] rhs - Value to invert.
 *
 * @return Result of bitwise NOT.
 */
template<typename E, EnableIfBitmaskOps<E> = 0>
constexpr E operator~(E rhs) noexcept
{
	using U = std::underlying_type_t<E>;

	rhs = static_cast<E>(~static_cast<U>(rhs));

	return rhs;
}

/**
 * Checks if all bits in rhs are set in lhs.
 *
 * @param[in] lhs - Value to test.
 * @param[in] rhs - Bitmask to check.
 *
 * @return true if all bits in rhs are set in lhs.
 */
template<typename E, EnableIfBitmaskOps<E> = 0>
[[nodiscard]] constexpr bool has_flag(E lhs, E rhs) noexcept
{
	return (lhs & rhs) == rhs;
}

/**
 * Sets or clears a flag in a bitmask.
 *
 * @param[in,out] flags - Bitmask to modify.
 * @param[in] flag - Flag to set or clear.
 * @param[in] set - true to set, false to clear.
 */
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

/*
################################################################################
Path Utilities
################################################################################
*/

/**
 * Extracts the file name from a path at compile-time.
 *
 * Returns a pointer into the input string. The input must be a string literal
 * (e.g., __FILE__) for the result to have valid lifetime. For owned strings,
 * use relative_path_storage().
 *
 * @param[in] path - Null-terminated path string (e.g., __FILE__).
 *
 * @return Pointer to the file name portion of the path.
 *
 * @warning Returns pointer into input string. Input must be a string literal
 *          for the returned pointer to remain valid.
 *
 * Example usage:
 * ```cpp
 * const char* filename = relative_path(__FILE__); // "CommonInclude.h"
 * ```
 */
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

/**
 * Extracts the file name from a path and stores it in a StackString.
 *
 * Unlike relative_path(), this returns an owned copy in a StackString, so the
 * input can be any null-terminated string (not just string literals).
 *
 * @param[in] path - Null-terminated path string.
 *
 * @return StackString containing the file name.
 */
[[nodiscard]] constexpr auto relative_path_storage(const char* path) noexcept
{
	StackString ret;

	ret.push_back(relative_path(path));

	return ret;
}

/*
################################################################################
String Utilities
################################################################################
*/

/**
 * Converts a string literal to a StackString at compile-time.
 *
 * Creates a StackString from a string literal. Useful for compile-time string
 * manipulation.
 *
 * @param[in] str - Null-terminated string literal.
 *
 * @return StackString containing the string.
 *
 * Example usage:
 * ```cpp
 * auto str = to_stack_string("hello world");
 * // str is StackString containing "hello world"
 * ```
 */
[[nodiscard]] constexpr auto to_stack_string(const char* str) noexcept
{
	StackString ret;

	ret.push_back(str);

	return ret;
}
