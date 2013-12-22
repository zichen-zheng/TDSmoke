#ifndef __GAUSS_FILTER_2D_H__
#define __GAUSS_FILTER_2D_H__

void GaussFilter2d(float *image, long width, long height,
    float sigma, int numsteps = 4);

#endif /* __GAUSS_FILTER_2D_H__ */
