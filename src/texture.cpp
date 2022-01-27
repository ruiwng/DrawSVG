#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

Color MipLevel::get_color_with_premuliply_alpha(size_t x, size_t y) {
  size_t pos = 4 * (x + y * width);
  float dst[4];
  uint8_to_float(dst, &texels[pos]);
  return Color(dst[0] * dst[3], dst[1] * dst[3], dst[2] * dst[3], dst[3]);
}

Color MipLevel::get_color(size_t x, size_t y) {
  size_t pos = 4 * (x + y * width);
  float dst[4];
  uint8_to_float(dst, &texels[pos]);
  return Color(dst[0], dst[1], dst[2], dst[3]);
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors (JUST AS A PLACEHOLDER)
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];
  /*
    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  */
    MipLevel& lastMipLevel = tex.mipmap[i - 1];
    for(size_t i = 0; i < mip.height; ++i) {
      for(size_t j = 0; j < mip.width; ++j) {
        Color color(0.0, 0.0, 0.0, 0.0);
        size_t x0 = j << 1;
        size_t y0 = i << 1;
        color += lastMipLevel.get_color(x0, y0);
        color += lastMipLevel.get_color(x0, y0 + 1);
        color += lastMipLevel.get_color(x0 + 1, y0);
        color += lastMipLevel.get_color(x0 + 1, y0 + 1);
        color *= 0.25;
        float_to_uint8(&mip.texels[4 * (j + i * mip.width)], &color.r);
      }
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 6: Implement nearest neighbour interpolation
  if (level >= tex.mipmap.size()) {
    return Color(1, 0, 1, 1);
  }
  // return Color(1.0, 0.0, 0.0, 1.0);
  // return magenta for invalid level
  MipLevel& mipLevel = tex.mipmap[level];
  size_t x = (mipLevel.width - 1) * u;
  size_t y = (mipLevel.height - 1) * v;
  return mipLevel.get_color_with_premuliply_alpha(x, y);
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering
  if (level >= tex.mipmap.size()) {
    return Color(1, 0, 1, 1);
  }

  MipLevel& mipLevel = tex.mipmap[level];
  float x = (mipLevel.width - 1) * u;
  float y = (mipLevel.height - 1) * v;
  size_t sx0 = floor(x);
  size_t sy0 = floor(y);
  size_t sx1 = sx0 + 1;
  if(sx1 >= mipLevel.width - 1) {
    sx1 = sx0;
  }
  size_t sy1 = sy0 + 1;
  if(sy1 >= mipLevel.height - 1) {
    sy1 = sy0;
  }
  float x_weight = x - sx0;
  float y_weight = y - sy0;

  // return magenta for invalid level
  auto color1 = mipLevel.get_color_with_premuliply_alpha(sx0, sy0);
  auto color2 = mipLevel.get_color_with_premuliply_alpha(sx0, sy1);
  auto color3 = mipLevel.get_color_with_premuliply_alpha(sx1, sy0);
  auto color4 = mipLevel.get_color_with_premuliply_alpha(sx1, sy1);

  return (color1 * (1.0 - y_weight) + color2 * y_weight) * (1.0 - x_weight) +  
  (color3 * (1.0 - y_weight) + color4 * y_weight) * x_weight;
}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering

  float numLevel = fmin(log2(1.0 / u_scale), log2(1.0 / v_scale));
  float level0 = floor(numLevel);
  float weight = numLevel - level0;
  size_t sLevel0 = fmax(0, fmin(level0, tex.mipmap.size() - 1));
  float level1 = ceil(numLevel);
  size_t sLevel1 = fmax(0, fmin(level1, tex.mipmap.size() - 1));

  Color color1 = sample_bilinear(tex, u, v, sLevel0);
  Color color2;
  if(sLevel0 != sLevel1) {
    color2 = sample_bilinear(tex, u, v, sLevel1);
  } else {
    color2 = color1;
  }
  return color1 * (1.0 - weight) + color2 * weight;
  // return magenta for invalid level
}

} // namespace CMU462
