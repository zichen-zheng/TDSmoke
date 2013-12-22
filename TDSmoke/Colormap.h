#ifndef __COLORMAP_H__
#define __COLORMAP_H__

#include <cmath>
#include <iostream>
#include <cassert>

#include "TDSmoke/MathUtilities.h"

class Colormap
{
public:

  // Colormaps supported by this class
  enum ColorScheme 
  {
    MATLAB_JET,
    //MATLAB_HSV,
    MATLAB_HOT,
    MATLAB_COOL,
    MATLAB_SPRING,
    MATLAB_SUMMER,
    MATLAB_AUTUMN,
    MATLAB_WINTER,
    MATLAB_GRAY,
    MATLAB_BONE,
    MATLAB_COPPER,
    MATLAB_PINK
  };  

  // Generate a color map with a given number of samples
  Colormap( const ColorScheme& color_scheme, const int& num_samples );
  
  // Deallocate any memory used by this object
  ~Colormap();

  // Set the number of samples in the color map
  void changeNumSamples( const int& num_samples );

  // Change the color map
  void changeColormap( const ColorScheme& color_scheme );

  // Change to the 'next' supported color map
  void incrementColormap();
  // Change to the 'previous' supported color map
  void decrementColormap();
  // TODO: refactor the above two methods so as not to repeat code.

  // Returns the number of samples of the current color map
  int getNumSamples() const;

  // Returns the current color map in use
  ColorScheme getColorScheme() const;

  // Given a scalar in [0,1.0], returns the corresponding color
  void getColorByDensity( const scalar& rho, scalar& r, scalar& g, scalar& b ) const;
  
  // Given an integer in [0,num_samples], returns the corresponding color
  void getColorByIndex( const int& i, scalar& r, scalar& g, scalar& b ) const;

  // Given a scalar in [0,1.0], returns the red channel
  scalar getRByDensity( const scalar& rho ) const;
  // Given a scalar in [0,1.0], returns the green channel
  scalar getGByDensity( const scalar& rho ) const;
  // Given a scalar in [0,1.0], returns the blue channel
  scalar getBByDensity( const scalar& rho ) const;  
  
  // Given an integer in [0,num_samples], returns the red channel
  scalar getRByIndex( const int& i ) const;
  // Given an integer in [0,num_samples], returns the green channel
  scalar getGByIndex( const int& i ) const;
  // Given an integer in [0,num_samples], returns the blue channel
  scalar getBByIndex( const int& i ) const;

  // Prints the values of the current color map
  friend std::ostream& operator<<( std::ostream& os, Colormap& rhs );

private:

  void generateColormap( const ColorScheme& color_scheme );
  void generateMatlabJet();
  void generateMatlabHot();
  void generateMatlabCool();
  void generateMatlabSpring();
  void generateMatlabSummer();
  void generateMatlabAutumn();
  void generateMatlabWinter();
  void generateMatlabGray();
  void generateMatlabBone();
  void generateMatlabCopper();
  void generateMatlabPink();

  int m_num_samples;
  ColorScheme m_color_scheme;

  // TODO: Move these to eigen types, clean up generation
  //  code accordingly.
  scalar* m_r;
  scalar* m_g;
  scalar* m_b;

};

#endif
