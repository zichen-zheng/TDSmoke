// YImage.h
//

#ifdef PNGOUT

#ifndef __YImage_h__
#define __YImage_h__

// file loading/saving automatically picks up changes to this struct.
// The possibilities are: ARGB, ABGR, RGBA, BGRA.

#include <iostream>
#include <cassert>
#include <cstdlib>

#include <png.h>

// for: #define offsetof(TYPE, MEMBER) ((size_t) &((TYPE *)0)->MEMBER)
// we use this to determine the pixel format on the fly
#include <cstddef>

struct YPixel
{
  unsigned char r;
  unsigned char g;
  unsigned char b;
  unsigned char a;
};

class YImage
{
public:
  YImage();

  YImage( const YImage& );

  virtual ~YImage();

  YImage& operator=( const YImage& );
  
  YPixel* data();
  
  const YPixel* data() const;
  
  YPixel& at( int i, int j );
  
  const YPixel& at( int i, int j ) const;

  int width() const;
  
  int height() const;  

  void resize(int width, int height);
  
  // Sets to greyscale by averaging pixels. TODO: Use the weighting that emphasizes green...
  void greyscale();
  
  // Flips the image vertically.
  void flip();
  
  // Flips the image horizontally.
  void mirror();

  bool save(const char* fname) const;

  bool load(const char* fname);


protected:
  int m_width;
  int m_height;
  YPixel* m_data;
};

#endif

#endif
