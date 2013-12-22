#include <math.h>
#include "GaussFilter2d.h"

/**
 * Fast 2D Gaussian convolution IIR approximation
 * @param image = the image data, modified in-place
 * @param width = width of image
 * @param height = heigth of image
 * @param sigma = the standard deviation of the Gaussian in pixels
 * @param numsteps = number of timesteps, more steps implies better accuracy
 *
 * Implements the fast Gaussian convolution algorithm of Alvarez and Mazorra,
 * where the Gaussian is approximated by a cascade of first-order infinite
 * impulsive response (IIR) filters.  Boundaries are handled with half-sample
 * symmetric extension.
 *
 * Gaussian convolution is approached as approximating the heat equation and
 * each timestep is performed with an efficient recursive computation.  Using
 * more steps yields a more accurate approximation of the Gaussian.  A
 * reasonable default value for \c numsteps is 4.
 *
 * The data is assumed to be ordered such that
 *   image[x + width*y] = pixel value at (x,y).
 *
 * Reference:
 * Alvarez, Mazorra, "Signal and Image Restoration using Shock Filters and
 * Anisotropic Diffusion," SIAM J. on Numerical Analysis, vol. 31, no. 2,
 * pp. 590-605, 1994.
 */

void GaussFilter2d(float *image, long width, long height,
        float sigma, int numsteps) {
    const long numpixels = width*height;
    double lambda, dnu;
    float nu, boundaryscale, postscale;
    float *ptr;
    long i, x, y;
    int step;
    
    if(sigma <= 0 || numsteps < 0)
        return;
    
    lambda = (sigma*sigma)/(2.0*numsteps);
    dnu = (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
    nu = (float)dnu;
    boundaryscale = (float)(1.0/(1.0 - dnu));
    postscale = (float)(pow(dnu/lambda,2*numsteps));
    
    /* Filter horizontally along each row */
    for(y = 0; y < height; y++) {
        for(step = 0; step < numsteps; step++) {
            ptr = image + width*y;
            ptr[0] *= boundaryscale;
            
            /* Filter rightwards */
            for(x = 1; x < width; x++)
                ptr[x] += nu*ptr[x - 1];
            
            ptr[x = width - 1] *= boundaryscale;
            
            /* Filter leftwards */
            for(; x > 0; x--)
                ptr[x - 1] += nu*ptr[x];
        }
    }
    
    /* Filter vertically along each column */
    for(x = 0; x < width; x++) {
        for(step = 0; step < numsteps; step++) {
            ptr = image + x;
            ptr[0] *= boundaryscale;
            
            /* Filter downwards */
            for(i = width; i < numpixels; i += width)
                ptr[i] += nu*ptr[i - width];
            
            ptr[i = numpixels - width] *= boundaryscale;
            
            /* Filter upwards */
            for(; i > 0; i -= width)
                ptr[i - width] += nu*ptr[i];
        }
    }
    
    for(i = 0; i < numpixels; i++)
        image[i] *= postscale;
    
    return;
}
