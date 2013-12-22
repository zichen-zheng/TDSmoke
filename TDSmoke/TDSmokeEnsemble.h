#ifndef __TD_SMOKE_ENSEMBLE_H__
#define __TD_SMOKE_ENSEMBLE_H__

#include "SimulationEnsemble.h"
#include "MathUtilities.h"
#include "RenderingUtilities.h"
#include "TDSmokeSim.h"
#include "Colormap.h"


class TDSmokeEnsemble : public SimulationEnsemble
{
public:
    
    TDSmokeEnsemble(int N);
    
    virtual ~TDSmokeEnsemble();
    
    virtual void stepSystem( const scalar& dt );
    
    /////////////////////////////////////////////////////////////////////////////
    // Mouse and keyboard input handlers
    
    virtual void keyboard( const unsigned char& key, const int& x, const int& y );
    
    virtual void special( const int& key, const int& x, const int& y );
    
    virtual void mouse( const int& button, const int& state, const int& x, const int& y );
    
    virtual void motion( const int& x, const int& y );
    
    /////////////////////////////////////////////////////////////////////////////
    // Functions for use by OpenGL and GLUT
    
    virtual void initializeOpenGL();
    
    virtual void reshape( const int& w, const int& h );
    
    virtual void display();
    virtual void display_cell();
    virtual void display_particle();
    virtual void display_surface();
    
    virtual unsigned int getGlutDisplayMode() const;
    
    virtual const int& getWindowWidth() const;
    
    virtual const int& getWindowHeight() const;
    
    virtual const renderingutils::Color& getBackgroundColor() const;

    virtual void initMarkerAndTarget();
    virtual void updateMarkerDensities();
    virtual void updateTargetDensities();

    virtual void setVelocityPattern( int velocity_pattern );
    virtual void setDiffusion( scalar diff );
    virtual void setViscosity( scalar visc );
    virtual void setSmoothing( scalar sigma );
    virtual void setDrivingForceCoeff( scalar vf );
    virtual void setAttenuation( scalar vd );
    virtual void setGathering( scalar vg );
    virtual void setGatheringEnabled( bool vg_enabled );
    
    virtual void setMarker( int i, int j );
    virtual void setTarget( int i, int j );
    
private:
    
    void addSphere( const int& row, const int& col, const int& R, bool is_target = false);
    ArrayXb generateSolidMap(int n);
    
    int m_N;
    ArrayXb m_solids;
    TDSmokeSim * m_fluid_sim;
    
    renderingutils::Color m_bgcolor;
    
    Colormap m_colormap;
    
    int m_window_width;
    int m_window_height;
    
    bool m_left_drag;
    bool m_right_drag;
    
    int m_last_row;
    int m_last_col;
    
    bool m_display_debugging;
    bool m_display_surface;
    
    Eigen::Array<std::vector<int>, Eigen::Dynamic, Eigen::Dynamic> m_particle_cell_map;
};

#endif
