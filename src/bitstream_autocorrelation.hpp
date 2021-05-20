#ifndef _DAISY_BISTREAM_AUTOCORRELATION_HPP_
#define _DAISY_BISTREAM_AUTOCORRELATION_HPP_

#include <cmath>
#include <cstdint>
#include <type_traits>
#include <vector>

// smallest power of 2 that fits n
template <typename T> constexpr T smallest_pow2(T n, T m = 1)
{
	return (m < n) ? smallest_pow2(n, m << 1) : m;
}

std::uint32_t count_bits(std::uint32_t i)
{
	// GCC only!!!
	return __builtin_popcount(i);
}

std::uint64_t count_bits(std::uint64_t i)
{
	// GCC only!!!
	return __builtin_popcountll(i);
}

template <typename T = std::uint32_t> struct bitstream
{
	static_assert(std::is_unsigned<T>::value, "T must be unsigned");
	static constexpr auto nbits = 8 * sizeof(T);

	bitstream(std::size_t size_)
	{
		size = smallest_pow2(size_);
		array_size = size / nbits;
		bits.resize(array_size, 0);
	}

	void clear() { std::fill(bits.begin(), bits.end(), 0); }

	void set(std::uint32_t i, bool val)
	{
		auto mask = 1 << (i % nbits);
		auto& ref = bits[i / nbits];
		ref ^= (-T(val) ^ ref) & mask;
	}

	bool get(std::uint32_t i) const
	{
		auto mask = 1 << (i % nbits);
		return (bits[i / nbits] & mask) != 0;
	}

	template <typename F> void auto_correlate(std::size_t start_pos, F f)
	{
		auto mid_array = (array_size / 2) - 1;
		auto mid_pos = size / 2;
		auto index = start_pos / nbits;
		auto shift = start_pos % nbits;

		for (auto pos = start_pos; pos != mid_pos; ++pos)
		{
			auto* p1 = bits.data();
			auto* p2 = bits.data() + index;
			auto count = 0;

			if (shift == 0)
			{
				for (auto i = 0; i != mid_array; ++i)
					count += count_bits(*p1++ ^ *p2++);
			}
			else
			{
				auto shift2 = nbits - shift;
				for (auto i = 0; i != mid_array; ++i)
				{
					auto v = *p2++ >> shift;
					v |= *p2 << shift2;
					count += count_bits(*p1++ ^ v);
				}
			}
			++shift;
			if (shift == nbits)
			{
				shift = 0;
				++index;
			}

			f(pos, count);
		}
	}

	std::vector<T> bits;
	std::size_t size;
	std::size_t array_size;
};

struct zero_cross
{
	bool operator()(float s)
	{
		if (s < -0.1f)
			y = 0;
		else if (s > 0.0f)
			y = 1;
		return y;
	}

	bool y = 0;
};

struct noise
{
	float operator()() const { return (float(rand()) / (RAND_MAX / 2)) - 1.0; }
};

constexpr double pi() { return std::atan(1) * 4; }
constexpr auto sps = 44100;
constexpr auto minFreq = 50.0;
constexpr auto maxFreq = 1200.0;

constexpr float minPeriod = float(sps) / maxFreq;
constexpr float maxPeriod = float(sps) / minFreq;

constexpr std::size_t bufferSize =
    smallest_pow2<std::size_t>(std::ceil(maxPeriod)) * 2;

class BitstreamCorrelation
{
private:
	float m_Buffer[bufferSize] = {0};
	std::size_t m_Index = 0;

public:
	void addSample(float value) { m_Buffer[m_Index++] = value; }

	/// Guesses the frequency whenever the buffer fills up. This can be greatly
	/// improved with a circular buffer, but I just wanted to get the basic code
	/// from
	/// https://github.com/cycfi/bitstream_autocorrelation/blob/master/bcf2.cpp
	/// working with minimal modification.
	float guessFrequency()
	{
		if (m_Index != bufferSize)
		{
			return 0.0f;
		}

		m_Index = 0;

		bitstream<> bin(bufferSize);

		zero_cross zc;

		for (auto i = 0; i != bufferSize; ++i)
		{
			bin.set(i, zc(m_Buffer[i]));
		}

		std::uint32_t maxCount = 0;
		std::uint32_t minCount = UINT32_MAX;
		std::size_t estIndex = 0;
		std::vector<std::uint32_t> corr(bufferSize / 2);
		bin.auto_correlate(
		    minPeriod,
		    [&corr, &maxCount, &minCount, &estIndex](auto pos, auto count) {
			    corr[pos] = count;
			    maxCount = std::max<std::uint32_t>(maxCount, count);
			    if (count < minCount)
			    {
				    minCount = count;
				    estIndex = pos;
			    }
		    });
		////////////////////////////////////////////////////////////////////////////
		// Handle harmonics
		auto subThreshold = 0.15 * maxCount;
		int maxDiv = estIndex / minPeriod;
		for (int div = maxDiv; div != 0; div--)
		{
			bool allStrong = true;
			float mul = 1.0f / div;

			for (int k = 1; k != div; k++)
			{
				int subPeriod = k * estIndex * mul;
				if (corr[subPeriod] > subThreshold)
				{
					allStrong = false;
					break;
				}
			}

			if (allStrong)
			{
				estIndex = estIndex * mul;
				break;
			}
		}

		////////////////////////////////////////////////////////////////////////////
		// Estimate the pitch

		// Get the start edge
		float prev = 0;
		auto startEdge = m_Buffer;
		for (; *startEdge <= 0.0f; ++startEdge)
		{
			prev = *startEdge;
		}
		auto dy = *startEdge - prev;
		auto dx1 = -prev / dy;

		// Get the next edge
		auto nextEdge = m_Buffer + estIndex - 1;
		for (; *nextEdge <= 0.0f; ++nextEdge)
		{
			prev = *nextEdge;
		}
		dy = *nextEdge - prev;
		auto dx2 = -prev / dy;

		float numSamples = (nextEdge - startEdge) + (dx2 - dx1);
		return sps / numSamples;
	}
};
#endif