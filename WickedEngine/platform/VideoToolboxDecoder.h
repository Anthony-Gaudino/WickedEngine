#pragma once

#include "wiPlatform.h"

#ifdef PLATFORM_MACOS

#include <cstdint>
#include <memory>

#include "wiVector.h"

namespace wi::video
{
	struct Video;

	namespace apple
	{
		struct Decoder;

		std::shared_ptr<Decoder> CreateDecoder(uint32_t width, uint32_t height, uint32_t padded_width, uint32_t padded_height, const wi::vector<uint8_t>& sps, const wi::vector<uint8_t>& pps);
		bool DecodeFrame(Decoder& decoder, const uint8_t* bitstream, size_t size, bool keyframe, wi::vector<uint8_t>& rgba_out, uint32_t display_width, uint32_t display_height);
		void Reset(Decoder& decoder);
	}
}

#endif // PLATFORM_MACOS
