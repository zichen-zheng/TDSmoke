#ifndef __TWO_DIMENSIONAL_DISPLAY_CONTROLLER_H__
#define __TWO_DIMENSIONAL_DISPLAY_CONTROLLER_H__

#include <cassert>
#include "RenderingUtilities.h"

// Sets orthographic projection for a 2D scene, handles mouse and keyboard
// input to zoom, track, etc. 
class TwoDimensionalDisplayController
{
public:
  TwoDimensionalDisplayController( const int& width, const int& height, const scalar& scalefactor = 1.0 );
  
  void reshape( const int& w, const int& h );
  
  void keyboard( const unsigned char& key, const int& x, const int& y );
  
  void special( const int& key, const int& x, const int& y );

  void mouse( const int& button, const int& state, const int& x, const int& y );
  
  void motion( const int& x, const int& y );

  void translateView( const scalar& dx, const scalar& dy );
  
  void zoomView( const scalar& dx, const scalar& dy );

  const int& getWindowWidth() const;  
  const int& getWindowHeight() const;

private:
  // Width of the window in pixels
  int m_window_width;
  // Height of the window in pixels
  int m_window_height;
  // Factor to 'zoom' in or out by
  scalar m_scale_factor;
  // Center of the display, x coord
  scalar m_center_x;
  // Center of the display, y coord
  scalar m_center_y;
  // True if the user is dragging the display left
  bool m_left_drag;
  // True if the user is dragging the display right
  bool m_right_drag;
  // Last position of the cursor in a drag, x coord
  int m_last_x;
  // Last position of the cursor in a drag, y coord
  int m_last_y;
};

#endif
