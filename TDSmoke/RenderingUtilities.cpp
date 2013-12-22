#include "RenderingUtilities.h"

namespace renderingutils
{

bool checkGLErrors()
{
  GLenum errCode;

  // TODO: Expand verbosity of error messages
  if( (errCode = glGetError()) != GL_NO_ERROR )
  {
    const GLubyte* errString = gluErrorString(errCode);
    std::cout << outputmod::startred << "OPENGL ERROR:" << outputmod::endred << std::flush;
    fprintf(stderr, " %s\n", errString);
    return false;
  }
  return true;
}    

void bressenhamLine( const int& inx0, const int& iny0, const int& inx1, const int& iny1, std::vector<int>& xpoints, std::vector<int>& ypoints )
{
  int x0 = inx0;
  int y0 = iny0;
  int x1 = inx1;
  int y1 = iny1;

  int Dx = x1 - x0; 
  int Dy = y1 - y0;
  int steep = (abs(Dy) >= abs(Dx));

  if (steep) 
  {
    std::swap(x0, y0);
    std::swap(x1, y1);
    // recompute Dx, Dy after swap
    Dx = x1 - x0;
    Dy = y1 - y0;
  }

  int xstep = 1;
  
  if (Dx < 0)
  {
    xstep = -1;
    Dx = -Dx;
  }
  
  int ystep = 1;
  
  if (Dy < 0) 
  {
    ystep = -1;		
    Dy = -Dy; 
  }
  
  const int& TwoDy = 2*Dy; 
  const int& TwoDyTwoDx = TwoDy - 2*Dx; // 2*Dy - 2*Dx
  int E = TwoDy - Dx; //2*Dy - Dx
  int y = y0;
  int xDraw, yDraw;	
  
  for (int x = x0; x != x1; x += xstep)
  {
    if (steep) 
    {
      xDraw = y;
      yDraw = x;
    } 
    else
    {			
      xDraw = x;
      yDraw = y;
    }
    
    // plot
    //plot(xDraw, yDraw);
    xpoints.push_back(xDraw);
    ypoints.push_back(yDraw);
    // next
    if (E > 0)
    {
      E += TwoDyTwoDx; //E += 2*Dy - 2*Dx;
      y = y + ystep;
    } 
    else 
    {
      E += TwoDy; //E += 2*Dy;
    }
  }
}  


///////////////////////////////////////////////////////////////////////////////
// COLOR METHODS

Color::Color()
: m_color(Vector3s::Zero())
{}

Color::Color( const scalar& r, const scalar& g, const scalar& b )
: m_color(r,g,b)
{
  assert( (m_color.array()>=0.0).all() );
  assert( (m_color.array()<=1.0).all() );
}

}
