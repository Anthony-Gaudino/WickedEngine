#include "platform/VideoToolboxDecoder.h"

#ifdef PLATFORM_MACOS

#import <CoreMedia/CoreMedia.h>
#import <CoreVideo/CoreVideo.h>
#import <VideoToolbox/VideoToolbox.h>

#include <algorithm>
#include <cstring>

namespace wi::video::apple
{
	struct Decoder
	{
		VTDecompressionSessionRef session = nullptr;
		CMVideoFormatDescriptionRef format = nullptr;
		uint32_t width = 0;
		uint32_t height = 0;
		uint32_t padded_width = 0;
		uint32_t padded_height = 0;
		CVPixelBufferRef latest = nullptr;
	};

	namespace
	{
		void ReleaseDecoderResources(Decoder& decoder)
		{
			if (decoder.latest)
			{
				CVPixelBufferRelease(decoder.latest);
				decoder.latest = nullptr;
			}
			if (decoder.session)
			{
				VTDecompressionSessionWaitForAsynchronousFrames(decoder.session);
				VTDecompressionSessionInvalidate(decoder.session);
				CFRelease(decoder.session);
				decoder.session = nullptr;
			}
			if (decoder.format)
			{
				CFRelease(decoder.format);
				decoder.format = nullptr;
			}
		}

		void DecoderDeleter(Decoder* decoder)
		{
			if (decoder == nullptr)
			{
				return;
			}
			ReleaseDecoderResources(*decoder);
			delete decoder;
		}

		void DecodeOutputCallback(
			void* decompressionOutputRefCon,
			void*,
			OSStatus status,
			VTDecodeInfoFlags,
			CVImageBufferRef imageBuffer,
			CMTime,
			CMTime)
		{
			Decoder* decoder = reinterpret_cast<Decoder*>(decompressionOutputRefCon);
			if (decoder == nullptr)
			{
				return;
			}
			if (decoder->latest)
			{
				CVPixelBufferRelease(decoder->latest);
				decoder->latest = nullptr;
			}
			if (status == noErr && imageBuffer != nullptr)
			{
				decoder->latest = CVPixelBufferRetain(imageBuffer);
			}
		}

		inline uint8_t ClampByte(int value)
		{
			return static_cast<uint8_t>(std::clamp(value, 0, 255));
		}

		bool ConvertNV12ToRGBA(CVPixelBufferRef pixelBuffer, wi::vector<uint8_t>& output, uint32_t display_width, uint32_t display_height)
		{
			if (pixelBuffer == nullptr)
			{
				return false;
			}

			const size_t required_size = size_t(display_width) * display_height * 4;
			output.resize(required_size);

			CVReturn lock_result = CVPixelBufferLockBaseAddress(pixelBuffer, kCVPixelBufferLock_ReadOnly);
			if (lock_result != kCVReturnSuccess)
			{
				return false;
			}

			const uint8_t* y_plane = static_cast<const uint8_t*>(CVPixelBufferGetBaseAddressOfPlane(pixelBuffer, 0));
			const uint8_t* uv_plane = static_cast<const uint8_t*>(CVPixelBufferGetBaseAddressOfPlane(pixelBuffer, 1));
			const size_t y_stride = CVPixelBufferGetBytesPerRowOfPlane(pixelBuffer, 0);
			const size_t uv_stride = CVPixelBufferGetBytesPerRowOfPlane(pixelBuffer, 1);

			for (uint32_t y = 0; y < display_height; ++y)
			{
				const uint8_t* row_y = y_plane + size_t(y) * y_stride;
				const uint8_t* row_uv = uv_plane + size_t(y / 2) * uv_stride;
				uint8_t* dst = output.data() + (size_t(y) * display_width * 4);
				for (uint32_t x = 0; x < display_width; ++x)
				{
					const uint8_t Y = row_y[x];
					const uint8_t U = row_uv[(x / 2u) * 2u + 0u];
					const uint8_t V = row_uv[(x / 2u) * 2u + 1u];

					const int C = int(Y) - 16;
					const int D = int(U) - 128;
					const int E = int(V) - 128;

					const int R = (298 * C + 409 * E + 128) >> 8;
					const int G = (298 * C - 100 * D - 208 * E + 128) >> 8;
					const int B = (298 * C + 516 * D + 128) >> 8;

					dst[x * 4 + 0] = ClampByte(R);
					dst[x * 4 + 1] = ClampByte(G);
					dst[x * 4 + 2] = ClampByte(B);
					dst[x * 4 + 3] = 255;
				}
			}

			CVPixelBufferUnlockBaseAddress(pixelBuffer, kCVPixelBufferLock_ReadOnly);
			return true;
		}
	}

	std::shared_ptr<Decoder> CreateDecoder(uint32_t width, uint32_t height, uint32_t padded_width, uint32_t padded_height, const wi::vector<uint8_t>& sps, const wi::vector<uint8_t>& pps)
	{
		if (sps.empty() || pps.empty())
		{
			return nullptr;
		}

		Decoder* raw_decoder = new Decoder();
		raw_decoder->width = width;
		raw_decoder->height = height;
		raw_decoder->padded_width = padded_width;
		raw_decoder->padded_height = padded_height;
		std::shared_ptr<Decoder> decoder(raw_decoder, DecoderDeleter);

		const uint8_t* parameter_sets[] = { sps.data(), pps.data() };
		size_t parameter_sizes[] = { sps.size(), pps.size() };
		CMVideoFormatDescriptionRef format = nullptr;
		OSStatus status = CMVideoFormatDescriptionCreateFromH264ParameterSets(
			kCFAllocatorDefault,
			2,
			parameter_sets,
			parameter_sizes,
			4,
			&format
		);
		if (status != noErr || format == nullptr)
		{
			return nullptr;
		}
		raw_decoder->format = format;

		CFMutableDictionaryRef destination_attributes = CFDictionaryCreateMutable(
			kCFAllocatorDefault,
			3,
			&kCFTypeDictionaryKeyCallBacks,
			&kCFTypeDictionaryValueCallBacks
		);

		if (destination_attributes == nullptr)
		{
			return nullptr;
		}

		int32_t pixel_format_value = static_cast<int32_t>(kCVPixelFormatType_420YpCbCr8BiPlanarVideoRange);
		int32_t buffer_width_value = static_cast<int32_t>(padded_width);
		int32_t buffer_height_value = static_cast<int32_t>(padded_height);
		CFNumberRef pixel_format = CFNumberCreate(kCFAllocatorDefault, kCFNumberSInt32Type, &pixel_format_value);
		CFNumberRef buffer_width = CFNumberCreate(kCFAllocatorDefault, kCFNumberSInt32Type, &buffer_width_value);
		CFNumberRef buffer_height = CFNumberCreate(kCFAllocatorDefault, kCFNumberSInt32Type, &buffer_height_value);

		CFDictionarySetValue(destination_attributes, kCVPixelBufferPixelFormatTypeKey, pixel_format);
		CFDictionarySetValue(destination_attributes, kCVPixelBufferWidthKey, buffer_width);
		CFDictionarySetValue(destination_attributes, kCVPixelBufferHeightKey, buffer_height);

		CFRelease(pixel_format);
		CFRelease(buffer_width);
		CFRelease(buffer_height);

		VTDecompressionOutputCallbackRecord callback = {};
		callback.decompressionOutputCallback = &DecodeOutputCallback;
		callback.decompressionOutputRefCon = raw_decoder;

		OSStatus session_status = VTDecompressionSessionCreate(
			kCFAllocatorDefault,
			raw_decoder->format,
			nullptr,
			destination_attributes,
			&callback,
			&raw_decoder->session
		);

		CFRelease(destination_attributes);

		if (session_status != noErr || raw_decoder->session == nullptr)
		{
			return nullptr;
		}

		VTSessionSetProperty(raw_decoder->session, kVTDecompressionPropertyKey_RealTime, kCFBooleanTrue);

		return decoder;
	}

	bool DecodeFrame(Decoder& decoder, const uint8_t* bitstream, size_t size, bool keyframe, wi::vector<uint8_t>& rgba_out, uint32_t display_width, uint32_t display_height)
	{
		if (decoder.session == nullptr || decoder.format == nullptr || bitstream == nullptr || size == 0)
		{
			return false;
		}

		if (decoder.latest)
		{
			CVPixelBufferRelease(decoder.latest);
			decoder.latest = nullptr;
		}

		CMBlockBufferRef block_buffer = nullptr;
		OSStatus status = CMBlockBufferCreateWithMemoryBlock(
			kCFAllocatorDefault,
			const_cast<uint8_t*>(bitstream),
			size,
			kCFAllocatorNull,
			nullptr,
			0,
			size,
			false,
			&block_buffer
		);
		if (status != noErr)
		{
			return false;
		}

		CMSampleBufferRef sample_buffer = nullptr;
		status = CMSampleBufferCreateReady(
			kCFAllocatorDefault,
			block_buffer,
			decoder.format,
			1,
			0,
			nullptr,
			1,
			reinterpret_cast<const size_t*>(&size),
			&sample_buffer
		);

		if (status != noErr || sample_buffer == nullptr)
		{
			if (block_buffer)
			{
				CFRelease(block_buffer);
			}
			return false;
		}

		CFArrayRef attachments = CMSampleBufferGetSampleAttachmentsArray(sample_buffer, true);
		if (attachments != nullptr && CFArrayGetCount(attachments) > 0)
		{
			CFMutableDictionaryRef attachment = (CFMutableDictionaryRef)CFArrayGetValueAtIndex(attachments, 0);
			if (attachment != nullptr)
			{
				if (keyframe)
				{
					CFDictionaryRemoveValue(attachment, kCMSampleAttachmentKey_NotSync);
				}
				else
				{
					CFDictionarySetValue(attachment, kCMSampleAttachmentKey_NotSync, kCFBooleanTrue);
				}
				CFDictionarySetValue(attachment, kCMSampleAttachmentKey_DisplayImmediately, kCFBooleanTrue);
			}
		}

		VTDecodeInfoFlags info_flags = 0;
		status = VTDecompressionSessionDecodeFrame(
			decoder.session,
			sample_buffer,
			0,
			nullptr,
			&info_flags
		);

		VTDecompressionSessionWaitForAsynchronousFrames(decoder.session);

		CFRelease(sample_buffer);
		CFRelease(block_buffer);

		if (status != noErr || decoder.latest == nullptr)
		{
			return false;
		}

		const bool convert_success = ConvertNV12ToRGBA(decoder.latest, rgba_out, std::min(display_width, decoder.width), std::min(display_height, decoder.height));
		CVPixelBufferRelease(decoder.latest);
		decoder.latest = nullptr;
		return convert_success;
	}

	void Reset(Decoder& decoder)
	{
		if (decoder.latest)
		{
			CVPixelBufferRelease(decoder.latest);
			decoder.latest = nullptr;
		}
		if (decoder.session)
		{
			VTDecompressionSessionWaitForAsynchronousFrames(decoder.session);
		}
	}
}

#endif // PLATFORM_MACOS
