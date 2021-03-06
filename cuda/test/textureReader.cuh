///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      since texture references are implicit static file-scope
//      variables, they cannot be passed as kernel parameters
//      (see http://forums.nvidia.com/lofiversion/index.php?t70630.html)
//      these structures can be used to work around this limitation
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef TEXTUREREADER_CUH
#define TEXTUREREADER_CUH


#include "complexTest.cuh"


namespace rpwa {

	namespace cuda {


		texture<float,  1, cudaReadModeElementType> floatTexture;
		struct floatTextureReader {

			typedef float texture_type;

			static __device__ float fetch(const int index) { return tex1Dfetch(floatTexture, index); }
                        
                        static HOST void bindTexture(const void*        deviceInData,
                                                     const unsigned int memSize)
                        { return checkCudaErrors(cudaBindTexture(0, floatTexture, deviceInData, memSize)); }

                        static HOST void unbindTexture() { return checkCudaErrors(cudaUnbindTexture(floatTexture)); }

                };


		texture<float2, 1, cudaReadModeElementType> float2Texture;
		struct float2TextureReader {

			typedef float2 texture_type;

			static __device__ float2 fetch(const int index) { return tex1Dfetch(float2Texture, index); }

                        static HOST void bindTexture(const void*        deviceInData,
                                                     const unsigned int memSize)
                        { return checkCudaErrors(cudaBindTexture(0, float2Texture, deviceInData, memSize)); }

                        static HOST void unbindTexture() { return checkCudaErrors(cudaUnbindTexture(float2Texture)); }

                };


		struct floatComplexTextureReader {

			typedef float2 texture_type;

			static __device__ complexTest<float2, float> fetch(const int index)
			{
				const float2 val = tex1Dfetch(float2Texture, index);
				return complexTest<float2, float>(val.x, val.y);
			}

                        static HOST void bindTexture(const void*        deviceInData,
                                                     const unsigned int memSize)
                        { return checkCudaErrors(cudaBindTexture(0, float2Texture, deviceInData, memSize)); }

                        static HOST void unbindTexture() { return checkCudaErrors(cudaUnbindTexture(float2Texture)); }

                };


		texture<float4, 1, cudaReadModeElementType> float4Texture;
		struct float4TextureReader {

			typedef float4 texture_type;

			static __device__ float4 fetch(const int index) { return tex1Dfetch(float4Texture, index); }

                        static HOST void bindTexture(const void*        deviceInData,
                                                     const unsigned int memSize)
                        { return checkCudaErrors(cudaBindTexture(0, float4Texture, deviceInData, memSize)); }
  
                        static HOST void unbindTexture() { return checkCudaErrors(cudaUnbindTexture(float4Texture)); }

                };


		texture<int4, 1, cudaReadModeElementType> int4Texture;
		struct doubleComplexTextureReader {

			typedef int4 texture_type;

			static __device__ complexTest<double2, double> fetch(const int index)
			{
				const int4 val = tex1Dfetch(int4Texture, index);
				return complexTest<double2, double>(__hiloint2double(val.y, val.x),
				                                    __hiloint2double(val.w, val.z));
			}

                        static HOST void bindTexture(const void*        deviceInData,
                                                     const unsigned int memSize)
                        { return checkCudaErrors(cudaBindTexture(0, int4Texture, deviceInData, memSize)); }
  
                        static HOST void unbindTexture() { return checkCudaErrors(cudaUnbindTexture(int4Texture)); }

                };


	}  // namespace cuda

}  // namespace rpwa


#endif  // TEXTUREREADER_CUH
