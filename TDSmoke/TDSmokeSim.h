#ifndef __TD_SMOKE_SIM_H__
#define __TD_SMOKE_SIM_H__

#include <iostream>

#include "MathUtilities.h"
#include "GaussFilter2d.h"

#define DEFAULT_DIFF    0.0001
#define DEFAULT_VISC    0.00000001
#define DEFAULT_SIGMA   1.0
#define DEFAULT_VF      8.0
#define DEFAULT_VD      1.1
#define DEFAULT_VG      1e-4
#define DEFAULT_VG_ENABLED 1

class TDSmokeSim
{
public:
    
    // Initialize the fluid simulation. Assumes rows == cols.
    //   rows: Number of rows in the Eulerian grid.
    //   cols: Number of cols in the Eulerian grid.
    //   diff: Diffusion coefficient for the passive markers.
    //   visc: Viscosity of the fluid. 
    TDSmokeSim( const int& rows, const int& cols, const scalar& diff = DEFAULT_DIFF, const scalar& visc = DEFAULT_VISC, const scalar& sigma = DEFAULT_SIGMA, const scalar& vf = DEFAULT_VF, const scalar& vd = DEFAULT_VD, const scalar& vg = DEFAULT_VG, const bool& vg_enabled = DEFAULT_VG_ENABLED);
    
    // Deallocates memory used during the simulation.
    ~TDSmokeSim();
    
    // Integrates the system forward in time by dt. 
    virtual void stepSystem( const scalar& dt );
    
    // Returns an array containing the marker densities. Note that the
    // boundary of this array is padded to handle boundary conditions. 
    const ArrayXs& getMarkerDensities() const;
    ArrayXs& getMarkerDensities();
    
    // Returns an array containing the target smoke densities. Note that the
    // boundary of this array is padded to handle boundary conditions. 
    const ArrayXs& getTargetDensities() const;
    ArrayXs& getTargetDensities();
    
    // Returns an array containing the horizontal components of the fluid velocity. 
    // Note that the boundary of this array is padded to handle boundary conditions. 
    ArrayXs& getHorizontalVelocities();
    const ArrayXs& getHorizontalVelocities() const;
    
    // Returns an array containing the vertical components of the fluid velocity. 
    // Note that the boundary of this array is padded to handle boundary conditions. 
    ArrayXs& getVerticalVelocities();
    const ArrayXs& getVerticalVelocities() const;
    
    // Returns the number of non-boundary rows.
    int physicalRows() const;
    // Returns the number of non-boundary columns.
    int physicalCols() const;
    
    // Updates blurred densities for smoke according to the current smoke densities
    void updateBlurredMarkerDensities();
    // Updates blurred target densities for smoke according to the current
    // target densities
    void updateBlurredTargetDensities();
    
    // Clears the density and velocity fields
    virtual void clear();
    
    // Prescribed velocity field for debugging
    virtual void setPrescribedVelocity(int p);

    // Diffusion coeff
    void setDiffusion(scalar diff);

    // Viscosity coeff
    void setViscosity(scalar visc);
    
    // Smoothing param
    void setSmoothing(scalar sigma);
    
    // Driving force coeff
    void setDrivingForceCoeff(scalar vf);
    
    // Attenuation coeff
    void setAttenuation(scalar vd);
    
    // Gathering coeff
    void setGathering(scalar vg);
    
    // Gathering enabled
    void setGatheringEnabled(bool vg_enabled);
    
    // Get/set verbose level
    bool verbose(void) { return VERBOSE; }
    void setVerbose(bool b) { VERBOSE = b; }

    void copyState( const TDSmokeSim& otherscene );

    virtual ArrayXb& getHasFluid();
    virtual const ArrayXb& getHasFluid() const;

protected:
    // Convenient grid access
    static scalar & d(ArrayXs * d, int i, int j) { return (*d)(i, j); }
    
    static scalar & u(ArrayXs * u, int i, int j) { return (*u)(i, j); }
    static scalar & v(ArrayXs * v, int i, int j) { return (*v)(i, j); }
    
    // Interpolator
    scalar interpolateD(ArrayXs * d, scalar i, scalar j);
    scalar interpolateU(ArrayXs * u, scalar i, scalar j);
    scalar interpolateV(ArrayXs * v, scalar i, scalar j);
    
    // Gaussian filter
    void gaussConv(int N, ArrayXs* d, ArrayXs* d0);
    
protected:
    // Time stepping utilities
    virtual void dens_step(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, ArrayXs * tg_d, ArrayXs * blr_tg_d, scalar diff, scalar vg, bool vg_enabled, scalar dt);
    virtual void vel_step(int N, ArrayXs* u, ArrayXs* v, ArrayXs* u0, ArrayXs* v0, ArrayXs* blr_d, ArrayXs* blr_tg_d, scalar visc, scalar vf, scalar vd, scalar dt);

    virtual void add_source(int N, ArrayXs * x, ArrayXs * x0, scalar dt);
    
    void applyForce(int N, ArrayXs* u, ArrayXs* v, ArrayXs* u0, ArrayXs* v0, ArrayXs* d, ArrayXs* tgd, scalar vf, scalar dt);

    virtual void diffuseD(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt);
    virtual void diffuseU(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt);
    virtual void diffuseV(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt);
    
    virtual void attenuate(int N, ArrayXs * x, ArrayXs * x0, scalar vt, scalar dt);
    
    virtual void advectD(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt);
    virtual void advectU(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt);
    virtual void advectV(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt);
    
    virtual void project(int N, ArrayXs * u, ArrayXs * v, ArrayXs * u0, ArrayXs * v0);
    
    virtual void gather(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * tg_d, ArrayXs * blr_tg_d, scalar vg, scalar dt);

    void SWAP(ArrayXs *& x1, ArrayXs *& x2);
    
protected:    
    // Width of the (square) grid
    int m_N;
    // Density of smoke particles: defined at center of cell
    ArrayXs m_d;
    // Density of target smoke particles
    ArrayXs m_tg_d;
    // Blurred smoke density
    ArrayXs m_blr_d;
    // Blurred target smoke density
    ArrayXs m_blr_tg_d;
    // Horizontal velocities: defined at left/right walls of cell
    ArrayXs m_u;
    // Vertical velocities:   defined at top/bottom walls of cell
    ArrayXs m_v;
    
    /*--------- Coefficients ---------*/
    // Diffusion coefficient for tracer particles
    scalar m_diff;
    // Viscosity of fluid
    scalar m_visc;
    // Force smoothing parameter (i.e. Gaussian filter parameter)
    scalar m_sigma;
    // Force coefficient: boost or weaken the driving force
    scalar m_vf;
    // Velocity attenuation coefficient: decreasing it results in 
    // smoke's overshooting the target
    scalar m_vd;
    // Gathering coefficient: increasing it reduces "stray" smoke
    scalar m_vg;
    // Whether gathering is enabled or not
    bool m_vg_enabled;
    /*------ End of Coefficients -----*/
    
    bool VERBOSE;

    ArrayXb m_all_ones;
};

#endif
