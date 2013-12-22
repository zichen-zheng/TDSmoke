#ifndef __RENDERING_UTILITIES_H__
#define __RENDERING_UTILITIES_H__

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <GL/glext.h>
#endif

#include <iostream>
#include <cstdio>

#include "MathDefines.h"
#include "StringUtilities.h"

#include <vector>

namespace renderingutils
{

// False => error
bool checkGLErrors();

// Code for Bressenham's adapted from wikipedia: http://en.wikipedia.org/wiki/Bresenham's_line_algorithm
void bressenhamLine( const int& inx0, const int& iny0, const int& inx1, const int& iny1, std::vector<int>& xpoints, std::vector<int>& ypoints );

// Class to represent a color
class Color
{
public:
  Color();
  Color( const scalar& r, const scalar& g, const scalar& b );

  inline
  const scalar& r() const
  {
    return m_color.x();
  }

  inline
  const scalar& g() const
  {
    return m_color.y();
  }

  inline
  const scalar& b() const
  {
    return m_color.z();
  }
  
private:
  Vector3s m_color;
};

}

#endif
