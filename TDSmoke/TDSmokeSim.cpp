#include "TDSmokeSim.h"
#include <Eigen/LU>

//#define VERBOSE (0)
#define MAX_ITERATIONS 20

TDSmokeSim::TDSmokeSim( const int& rows, const int& cols, const scalar& diff, const scalar& visc, const scalar& sigma, const scalar& vf, const scalar& vd, const scalar& vg, const bool& vg_enabled)
: m_N(rows)
, m_d(m_N + 2, m_N + 2)
, m_tg_d(m_N + 2, m_N + 2)
, m_blr_d(m_N + 2, m_N + 2)
, m_blr_tg_d(m_N + 2, m_N + 2)
, m_u(m_N + 2, m_N + 1)
, m_v(m_N + 1, m_N + 2)
, m_diff(diff)
, m_visc(visc)
, m_sigma(sigma)
, m_vf(vf)
, m_vd(vd)
, m_vg(vg)
, m_vg_enabled(vg_enabled)
, VERBOSE(false)
, m_all_ones(m_N + 2, m_N + 2) {
    assert(rows==cols);
    
    clear();
}

TDSmokeSim::~TDSmokeSim() {
    // destructor: do nothing
}

////////////////////////////////////////////////////////////////////////////////////
void TDSmokeSim::SWAP(ArrayXs *& x1, ArrayXs *& x2) {
    ArrayXs * t = x1;
    x1 = x2;
    x2 = t;
}

void TDSmokeSim::add_source(int N, ArrayXs * x, ArrayXs * x0, scalar dt) {
    // nothing to do. adding marker/velocity is handled in TDSmokeEnsemble.
    *x = *x0;
}

void TDSmokeSim::applyForce(int N, ArrayXs* u, ArrayXs* v, ArrayXs* u0, ArrayXs* v0, ArrayXs* d, ArrayXs* tgd, scalar vf, scalar dt) {
    assert(u->cols() == u0->cols());
    assert(u->rows() == u0->rows());
    assert(v->cols() == v0->cols());
    assert(v->rows() == v0->rows());

    scalar dx = 1.0 / N;
    u->setZero();
    v->setZero();
    
    for (int i = 0; i < N+2; i++) {
        for (int j = 0; j < N+1; j++) {
            // compute driving force on u direction
            scalar d_u = interpolateD(d, i, j+0.5);
            scalar tgd_u = interpolateD(tgd, i, j+0.5);
            scalar d_tgd_u = ((*tgd)(i,j+1) - (*tgd)(i,j)) / dx;
            scalar f_u = d_u * d_tgd_u / tgd_u;
            
            // compute driving force on v direction
            scalar d_v = interpolateD(d, j+0.5, i);
            scalar tgd_v = interpolateD(tgd, j+0.5, i);
            scalar d_tgd_v = ((*tgd)(j+1,i) - (*tgd)(j,i)) / dx;
            scalar f_v = d_v * d_tgd_v / tgd_v;
            
            // update velocities
            (*u)(i,j) = (*u0)(i,j) + dt * vf * f_u;
            (*v)(j,i) = (*v0)(j,i) + dt * vf * f_v;
        }
    }
}

void TDSmokeSim::diffuseD(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt) {
    scalar dx = 1.0 / N;
    // initial guess of solution
    (*x) = (*x0);
    // solve by Gauss-Seidel iterations
    int currIteration = 0;
    while (currIteration < MAX_ITERATIONS) {
        for (int i = 0; i < N+2; i++) {
            for (int j = 0; j < N+2; j++) {
                /*
                int x_factor = 4;
                scalar all_x = (*x)(i-1,j) + (*x)(i+1,j) + (*x)(i,j-1) + (*x)(i,j+1);
                (*x)(i,j) = (*x0)(i,j)*pow(dx,2.0) + diff*dt*all_x;
                (*x)(i,j) /= ( pow(dx,2.0) + x_factor*diff*dt );
                */
                scalar all_x = 0;
                int x_factor = 0;
                // solid boundary conditions
                if (i != 0) {
                    all_x += (*x)(i-1,j);
                    x_factor ++;
                }
                if (i != N+1) {
                    all_x += (*x)(i+1,j);
                    x_factor ++;
                }
                if (j != 0) {
                    all_x += (*x)(i,j-1);
                    x_factor ++;
                }
                if (j != N+1) {
                    all_x += (*x)(i,j+1);
                    x_factor ++;
                }
                // update solution
                (*x)(i,j) = (*x0)(i,j)*pow(dx,2.0) + diff*dt*all_x;
                (*x)(i,j) /= ( pow(dx,2.0) + x_factor*diff*dt );
            }
        }
        currIteration ++;
    }
}

void TDSmokeSim::diffuseU(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt) {
    scalar dx = 1.0 / N;
    // initial guess of solution
    (*x) = (*x0);
    // solve by Gauss-Seidel iterations
    int currIteration = 0;
    while (currIteration < MAX_ITERATIONS) {
        for (int i = 1; i < N+1; i++) {
            for (int j = 0; j < N+1; j++) {
                scalar all_x = (*x)(i-1,j) + (*x)(i+1,j);
                int x_factor = 2;
                // solid boundary conditions
                if (j != 0) {
                    all_x += (*x)(i,j-1);
                    x_factor ++;
                }
                if (j != N) {
                    all_x += (*x)(i,j+1);
                    x_factor ++;
                }
                // update solution
                (*x)(i,j) = (*x0)(i,j)*pow(dx,2.0) + diff*dt*all_x;
                (*x)(i,j) /= ( pow(dx,2.0) + x_factor*diff*dt );
            }
        }
        currIteration ++;
    }
    
}

void TDSmokeSim::diffuseV(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt) {
    scalar dx = 1.0 / N;
    // initial guess of solution
    (*x) = (*x0);
    // solve by Gauss-Seidel iterations
    int currIteration = 0;
    while (currIteration < MAX_ITERATIONS) {
        for (int i = 0; i < N+1; i++) {
            for (int j = 1; j < N+1; j++) {
                scalar all_x = (*x)(i,j-1) + (*x)(i,j+1);
                int x_factor = 2;
                // solid boundary conditions
                if (i != 0) {
                    all_x += (*x)(i-1,j);
                    x_factor ++;
                }
                if (i != N) {
                    all_x += (*x)(i+1,j);
                    x_factor ++;
                }
                // update solution
                (*x)(i,j) = (*x0)(i,j)*pow(dx,2.0) + diff*dt*all_x;
                (*x)(i,j) /= ( pow(dx,2.0) + x_factor*diff*dt );
            }
        }
        currIteration ++;
    }
}

void TDSmokeSim::attenuate(int N, ArrayXs * x, ArrayXs * x0, scalar vd, scalar dt) {
    // initial guess of solution
    (*x) = (*x0);
    // solve by Gauss-Seidel iterations
    int currIteration = 0;
    while (currIteration < MAX_ITERATIONS) {
        for (int i = 0; i < x0->rows(); i++) {
            for (int j = 0; j < x0->cols(); j++) {
                // update solution
                (*x)(i,j) /= (dt * vd + 1);
            }
        }
        currIteration ++;
    }
}

void TDSmokeSim::advectD(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt) {
    scalar h = 1.0 / N;
    for (int i = 1; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            // current (interpolated) velocity
            scalar currU = interpolateU(u, i, j-0.5);
            scalar currV = interpolateV(v, i-0.5, j);
            // time integration backward in time
            scalar dy = i*h - currV*dt;
            scalar dx = j*h - currU*dt;
            // interpolate density
            (*x)(i,j) = interpolateD(x0, dy*N, dx*N);
        }
    }
}

void TDSmokeSim::advectU(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt) {
    scalar h = 1.0 / N;
    for (int i = 1; i < N+1; i++) {
        for (int j = 0; j < N+1; j++) {
            // current (interpolated) velocity
            scalar currU = interpolateU(u, i, j);
            scalar currV = interpolateV(v, i-0.5, j+0.5);
            // time integration backward in time
            scalar dy = i*h - currV*dt;
            scalar dx = j*h - currU*dt;
            // interpolate velocity in u direction
            (*x)(i,j) = interpolateU(x0, dy*N, dx*N);
        }
    }
}

void TDSmokeSim::advectV(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt) {
    scalar h = 1.0 / N;
    for (int i = 0; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            // current (interpolated) velocity
            scalar currU = interpolateU(u, i+0.5, j-0.5);
            scalar currV = interpolateV(v, i, j);
            // time integration backward in time
            scalar dy = i*h - currV*dt;
            scalar dx = j*h - currU*dt;
            // interpolate velocity in v direction
            (*x)(i,j) = interpolateV(x0, dy*N, dx*N);
        }
    }
}

void TDSmokeSim::gather(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * tg_d, ArrayXs * blr_tg_d, scalar vg, scalar dt) {
    scalar dx = 1.0 / N;
    // initial guess of solution
    (*x) = (*x0);
    // solve by Gauss-Seidel iterations
    int currIteration = 0;
    while (currIteration < MAX_ITERATIONS) {
        for (int i = 0; i < N+2; i++) {
            for (int j = 0; j < N+2; j++) {
                scalar all_rhs_x = 0;
                scalar all_lhs_x = 0;
                // solid boundary conditions
                if (i != 0) {
                    all_rhs_x += ( (*x)(i-1,j)-(*tg_d)(i-1,j)+(*tg_d)(i,j) ) * (*blr_tg_d)(i,j);
                    all_lhs_x += (*blr_tg_d)(i,j);
                }
                if (i != N+1) {
                    all_rhs_x += ( (*x)(i+1,j)-(*tg_d)(i+1,j)+(*tg_d)(i,j) ) * (*blr_tg_d)(i+1,j) * (*x)(i+1,j);
                    all_lhs_x += (*blr_tg_d)(i+1,j) * (*x)(i+1,j);
                }
                if (j != 0) {
                    all_rhs_x += ( (*x)(i,j-1)-(*tg_d)(i,j-1)+(*tg_d)(i,j) ) * (*blr_tg_d)(i,j);
                    all_lhs_x += (*blr_tg_d)(i,j);
                }
                if (j != N+1) {
                    all_rhs_x += ( (*x)(i,j+1)-(*tg_d)(i,j+1)+(*tg_d)(i,j) ) * (*blr_tg_d)(i,j+1) * (*x)(i,j+1);
                    all_lhs_x += (*blr_tg_d)(i,j+1) * (*x)(i,j+1);
                }
                // update solution
                (*x)(i,j) = (*x0)(i,j)*pow(dx,2.0) + vg * dt * all_rhs_x;
                (*x)(i,j) /= ( pow(dx,2.0) + vg*dt*all_lhs_x );
            }
        }
        currIteration ++;
    }
}

void TDSmokeSim::project(int N, ArrayXs * u, ArrayXs * v, ArrayXs * u0, ArrayXs * v0) {
    // Step 1: solid boundary condition for velocity
    for (int i = 0; i < N+1; i++) {
        // on a solid wall
        if (i != 0) {
            (*u0)(i,0) = 0;
            (*u0)(i,N) = 0;
            (*v0)(0,i) = 0;
            (*v0)(N,i) = 0;
        }
        // inside solids
        (*u0)(0,i) = 0;
        (*u0)(N+1,i) = 0;
        (*v0)(i,0) = 0;
        (*v0)(i,N+1) = 0;
    }
    
    // Step 2: compute the divergence of the velocity field
    ArrayXs gradU = ArrayXs::Zero(N+2,N+2);
    for (int i = 1; i < N+1; i++) {
        for (int j = 1; j < N+1; j++) {
            scalar tmp_gradU_u = (*u0)(i,j) - (*u0)(i,j-1);
            scalar tmp_gradU_v = (*v0)(i,j) - (*v0)(i-1,j);
            gradU(i,j) = tmp_gradU_u + tmp_gradU_v;
        }
    }
    gradU *= N;
    
    // Step 3: compute the pressure field
    int currIteration = 0;
    // initial guess to the solution
    ArrayXs p = ArrayXs::Zero(N+2,N+2);
    scalar dx = 1.0 / N;
    // solve by Gauss-Seidel iterations
    while (currIteration < MAX_ITERATIONS) {
        for (int i = 1; i < N+1; i++) {
            for (int j = 1; j < N+1; j++) {
                scalar all_p = 0;
                int p_factor = 0;
                // solid boundary conditions for pressure
                if (i != 1) {
                    all_p += p(i-1,j);
                    p_factor ++;
                }
                if (j != 1) {
                    all_p += p(i,j-1);
                    p_factor ++;
                }
                if (i != N) {
                    all_p += p(i+1,j);
                    p_factor ++;
                }
                if (j != N) {
                    all_p += p(i,j+1);
                    p_factor ++;
                }
                // update pressure by solving the linear equation
                p(i,j) = ( all_p - gradU(i,j)*pow(dx,2.0) ) / p_factor;
            }
        }
        currIteration ++;
    }
    
    // Step 4: update the velocity field
    // for u
    (*u) = ArrayXs::Zero(N+2,N+1);  // initialize
    assert( u->rows()==u0->rows() && u->cols()==u0->cols() );
    for (int i = 1; i < N+1; i++)
        for (int j = 1; j < N; j++)
            (*u)(i,j) = (*u0)(i,j) - (p(i,j+1) - p(i,j)) / dx;
    // for v
    (*v) = ArrayXs::Zero(N+1,N+2);  // initialize
    assert( v->rows()==v0->rows() && v->cols()==v0->cols() );
    for (int i = 1; i < N; i++)
        for (int j = 1; j < N+1; j++)
            (*v)(i,j) = (*v0)(i,j) - (p(i+1,j) - p(i,j)) / dx;
}

scalar TDSmokeSim::interpolateD(ArrayXs * d, scalar i, scalar j) {
    assert(d->rows() == d->cols());
    
    int N = d->rows() - 2;

    if (i < 0.5) i = 0.5;
    else if (i > N+0.5) i = N+0.5;
    if (j < 0.5) j = 0.5;
    else if (j > N+0.5) j = N+0.5;
    
    int ii = (int) i;
    int jj = (int) j;

    // weights for bilinear interpolation
    scalar iFrac = 1.0 - ( i - (scalar)ii );
    scalar jFrac = 1.0 - ( j - (scalar)jj );

    return iFrac*jFrac*(*d)(ii,jj) + iFrac*(1.0-jFrac)*(*d)(ii,jj+1) + (1.0-iFrac)*jFrac*(*d)(ii+1,jj) + (1.0-iFrac)*(1.0-jFrac)*(*d)(ii+1,jj+1);
}

scalar TDSmokeSim::interpolateU(ArrayXs * u, scalar i, scalar j) {
    assert(u->rows() == u->cols()+1);

    int N = u->rows() - 2;

    if (i < 0.5) i = 0.5;
    else if (i > (scalar)N+0.5) i = (scalar)N+0.5;
    if (j < 0) j = 0;
    else if (j > N) j = N;
    
    int ii = (int) i;
    int jj = (int) j;
    if (jj == N) jj = N-1;

    // weights for bilinear interpolation
    scalar iFrac = 1.0 - ( i - (scalar)ii );
    scalar jFrac = 1.0 - ( j - (scalar)jj );

    return iFrac*jFrac*(*u)(ii,jj) + iFrac*(1.0-jFrac)*(*u)(ii,jj+1) + (1.0-iFrac)*jFrac*(*u)(ii+1,jj) + (1.0-iFrac)*(1.0-jFrac)*(*u)(ii+1,jj+1);
}

scalar TDSmokeSim::interpolateV(ArrayXs * v, scalar i, scalar j) {
    assert(v->rows()+1 == v->cols());

    int N = v->cols() - 2;

    if (i < 0) i = 0;
    else if (i > N) i = N;
    if (j < 0.5) j = 0.5;
    else if (j > N+0.5) j = N+0.5;
    
    int ii = (int) i;
    int jj = (int) j;
    if (ii == N) ii = N-1;

    // weights for bilinear interpolation
    scalar iFrac = 1.0 - ( i - (scalar)ii );
    scalar jFrac = 1.0 - ( j - (scalar)jj );

    return iFrac*jFrac*(*v)(ii,jj) + iFrac*(1.0-jFrac)*(*v)(ii,jj+1) + (1.0-iFrac)*jFrac*(*v)(ii+1,jj) + (1.0-iFrac)*(1.0-jFrac)*(*v)(ii+1,jj+1);
}

void TDSmokeSim::dens_step(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, ArrayXs * tg_d, ArrayXs * blr_tg_d, scalar diff, scalar vg, bool vg_enabled, scalar dt) {
    ArrayXs * outu = u;
    ArrayXs * outv = v;
    
    add_source(N, x, x0, dt); 
    
    SWAP(x0, x); 
    if (vg_enabled)
        gather(N, x, x0, tg_d, blr_tg_d, vg, dt); 
    else
        diffuseD(N, x, x0, diff, dt);
    
    SWAP(x0, x); 
    advectD(N, x, x0, u, v, dt);
    
    if (outu != u)
        *outu = *u;    
    if (outv != v)
        *outv = *v;
}

void TDSmokeSim::vel_step(int N, ArrayXs* u, ArrayXs* v, ArrayXs* u0, ArrayXs* v0, ArrayXs* blr_d, ArrayXs* blr_tg_d, scalar visc, scalar vf, scalar vd, scalar dt) {
    ArrayXs * outu = u;
    ArrayXs * outv = v;
    
    add_source(N, u, u0, dt); 
    add_source(N, v, v0, dt); 
    
    SWAP(u0, u); 
    SWAP(v0, v); 
    applyForce(N, u, v, u0, v0, blr_d, blr_tg_d, vf, dt);
    
    SWAP(u0, u); 
    SWAP(v0, v); 
    diffuseU(N, u, u0, visc, dt); 
    diffuseV(N, v, v0, visc, dt); 

    SWAP(u0, u);
    SWAP(v0, v);
    project(N, u, v, u0, v0);
    
    SWAP(u0, u); 
    SWAP(v0, v); 
    attenuate(N, u, u0, vd, dt); 
    attenuate(N, v, v0, vd, dt); 
    
    SWAP(u0, u); 
    SWAP(v0, v); 
    advectU(N, u, u0, u0, v0, dt); 
    advectV(N, v, v0, u0, v0, dt); 
    
    SWAP(u0, u);
    SWAP(v0, v);
    project(N, u, v, u0, v0);
    
    if (outu != u)
        *outu = *u;    
    if (outv != v)
        *outv = *v;
    
}

void TDSmokeSim::stepSystem( const scalar& dt) {
    if (VERBOSE) std::cout << "step" << std::endl;
    ArrayXs new_d(m_N + 2, m_N + 2);
    ArrayXs new_blr_d(m_N + 2, m_N + 2);
    ArrayXs new_u(m_N + 2, m_N + 1);
    ArrayXs new_v(m_N + 1, m_N + 2);
    
    new_d.setZero();
    new_blr_d.setZero();
    new_u.setZero();
    new_v.setZero();
    
    vel_step(m_N, &new_u, &new_v, &m_u, &m_v, &m_blr_d, &m_blr_tg_d, m_visc, m_vf, m_vd, dt);
    dens_step(m_N, &new_d, &m_d, &new_u, &new_v, &m_tg_d, &m_blr_tg_d, m_diff, m_vg, m_vg_enabled, dt);
    
    gaussConv(m_N, &new_blr_d, &m_d);
    
    m_d = new_d;
    m_blr_d = new_blr_d;
    m_u = new_u;
    m_v = new_v;
}

void TDSmokeSim::gaussConv(int N, ArrayXs* d, ArrayXs* d0) {
    assert(d0->rows() == N+2);
    assert(d0->cols() == N+2);
    assert(d->rows() == N+2);
    assert(d->cols() == N+2);
    
    float d_arr[(N+2)*(N+2)];
    
    // read density from ArrayXs to float*
    for (int i = 0; i < N+2; i++) {
        for (int j = 0; j < N+2; j++) {
            d_arr[j+i*(N+2)] = (*d0)(i,j);
        }
    }
    
    // pass density to Gaussian filter
    GaussFilter2d(d_arr, (long) N+2, (long) N+2, m_sigma);
    
    // read density from float* back to ArrayXs
    for (int i = 0; i < N+2; i++) {
        for (int j = 0; j < N+2; j++) {
            (*d)(i,j) = d_arr[j+i*(N+2)];
            
            // To ensure the density gradient is NOT 0 everywhere,
            // i.e. to avoid situations where densities may be constant
            // in some regions (in our case, densities are 0 in some
            // regions), we add some noise when densities are too close
            // to 0.
            if ((*d)(i,j) < 1e-60)
                (*d)(i,j) = (scalar) rand()/RAND_MAX * 1e-59;
            
            //std::cout << (*d)(i,j) << std::endl;
        }
    }
}

void TDSmokeSim::updateBlurredMarkerDensities() {
    ArrayXs new_blr_d(m_N + 2, m_N + 2);
    gaussConv(m_N, &new_blr_d, &m_d);
    m_blr_d = new_blr_d;
}

void TDSmokeSim::updateBlurredTargetDensities() {
    ArrayXs new_blr_tg_d(m_N + 2, m_N + 2);
    gaussConv(m_N, &new_blr_tg_d, &m_tg_d);
    m_blr_tg_d = new_blr_tg_d;
}

const ArrayXs& TDSmokeSim::getMarkerDensities() const { return m_d; }
ArrayXs& TDSmokeSim::getMarkerDensities() { return m_d; }

const ArrayXs& TDSmokeSim::getTargetDensities() const { return m_tg_d; }
ArrayXs& TDSmokeSim::getTargetDensities() { return m_tg_d; }

ArrayXs& TDSmokeSim::getHorizontalVelocities() { return m_u; }
const ArrayXs& TDSmokeSim::getHorizontalVelocities() const { return m_u;}

ArrayXs& TDSmokeSim::getVerticalVelocities() { return m_v; }
const ArrayXs& TDSmokeSim::getVerticalVelocities() const { return m_v; }

int TDSmokeSim::physicalRows() const { return m_N; }
int TDSmokeSim::physicalCols() const { return m_N; }

void TDSmokeSim::clear() {
    m_d.setZero();
    m_tg_d.setZero();
    m_blr_d.setZero();
    m_blr_tg_d.setZero();
    m_u.setZero();
    m_v.setZero();
    m_all_ones.setOnes();
}

void TDSmokeSim::setPrescribedVelocity(int p) {
    switch (p) {
        case 0:
            break;
        case 1:
            for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
                for (int j = m_N * 0.2; j <= m_N * 0.8; j++)
                    m_u(i, j) = m_v(i, j) = 0.8;
            break;
        case 2:
            for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
                m_u(i, i) = m_v(i, i) = 0.8;
            break;
        case 3:
            for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
                m_v(i, m_N * 0.2 + 1) = -1.6,
                m_v(i, m_N * 0.8 + 1) = -1.6;
            break;
        case 4:
            for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
                m_u(m_N, i) = 2.4;
            break;
        case 5:
            for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
                for (int j = m_N * 0.83; j <= m_N; j++)
                    m_u(j, i) = 2.4;
            break;
        case 8:
            for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
                m_v(i, m_N / 2) = 0.8;
            break;
        case 9:
            for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
                m_v(i, m_N / 2) = -0.8;
            break;
    }
}


void TDSmokeSim::copyState( const TDSmokeSim& otherscene ) {
    m_diff = otherscene.m_diff;
    m_visc = otherscene.m_visc;
    m_N = otherscene.m_N;
    m_d = otherscene.m_d;
    m_u = otherscene.m_u;
    m_v = otherscene.m_v;
}

ArrayXb& TDSmokeSim::getHasFluid()  {return m_all_ones; }
const ArrayXb& TDSmokeSim::getHasFluid() const { return m_all_ones; }

void TDSmokeSim::setDiffusion(scalar diff) { m_diff = diff; }
void TDSmokeSim::setViscosity(scalar visc) { m_visc = visc; }
void TDSmokeSim::setSmoothing(scalar sigma) { m_sigma = sigma; }
void TDSmokeSim::setDrivingForceCoeff(scalar vf) { m_vf = vf; }
void TDSmokeSim::setAttenuation(scalar vd) { m_vd = vd; }
void TDSmokeSim::setGathering(scalar vg) { m_vg = vg; }
void TDSmokeSim::setGatheringEnabled(bool vg_enabled) { m_vg_enabled = vg_enabled; }


