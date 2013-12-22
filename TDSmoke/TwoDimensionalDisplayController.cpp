#include "TwoDimensionalDisplayController.h"

TwoDimensionalDisplayController::TwoDimensionalDisplayController( const int& width, const int& height, const scalar& scalefactor )
: m_window_width(width)
, m_window_height(height)
, m_scale_factor(scalefactor)
, m_center_x(0.0)
, m_center_y(0.0)
, m_left_drag(false)
, m_right_drag(false)
, m_last_x(0)
, m_last_y(0)
{}

void TwoDimensionalDisplayController::reshape( const int& w, const int& h )
{
    assert( renderingutils::checkGLErrors() );
    // Record the new width and height
    m_window_width = w;
    m_window_height = h;
    // Reset the coordinate system before modifying
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // Set the coordinate system to achieve the desired zoom level, center
    const scalar& ratio = ((scalar)h)/((scalar)w);
    gluOrtho2D(m_center_x-m_scale_factor/ratio,m_center_x+m_scale_factor/ratio,m_center_y-m_scale_factor,m_center_y+m_scale_factor);
    // Set the viewport to be the entire window
    glViewport(0,0,w,h);
    // Render the scene
    glutPostRedisplay();
    assert( renderingutils::checkGLErrors() );
}

void TwoDimensionalDisplayController::keyboard( const unsigned char& key, const int& x, const int& y )
{
    // Zoom out
    if( key == '-' || key == '_' )
    {
        m_scale_factor += 0.1;
        reshape(m_window_width,m_window_height);
    }
    // Zoom in
    else if( key == '=' || key == '+' )
    {
        m_scale_factor = std::max(0.1,m_scale_factor-0.1);
        reshape(m_window_width,m_window_height);
    }
}

void TwoDimensionalDisplayController::special( const int& key, const int& x, const int& y )
{
    // Translate the viewport up
    if( GLUT_KEY_UP == key )
    {
        m_center_y += 0.1;
        reshape(m_window_width,m_window_height);
    }
    // Translate the viewport down
    else if( GLUT_KEY_DOWN == key )
    {
        m_center_y -= 0.1;
        reshape(m_window_width,m_window_height);
    }
    // Translate the viewport left
    else if( GLUT_KEY_LEFT == key )
    {
        m_center_x -= 0.1;
        reshape(m_window_width,m_window_height);
    }
    // Translate the viewport right
    else if( GLUT_KEY_RIGHT == key )
    {
        m_center_x += 0.1;
        reshape(m_window_width,m_window_height);
    }
}

void TwoDimensionalDisplayController::mouse( const int& button, const int& state, const int& x, const int& y )
{
    // Record that the left mouse button was pressed (begin left drag)
    if( !m_right_drag && button == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
    {
        m_left_drag = true;
        m_last_x = x;
        m_last_y = y;
    }
    // Record that the left mouse button was released (end left drag)
    if( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
    {
        m_left_drag = false;
    }
    
    // Record that the right mouse button was pressed (begin right drag)
    if( !m_left_drag && button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN )
    {
        m_right_drag = true;
        m_last_x = x;
        m_last_y = y;
    }
    // Record that the right mouse button was released (end right drag)
    if( button == GLUT_RIGHT_BUTTON && state == GLUT_UP )
    {
        m_right_drag = false;
    }
}

void TwoDimensionalDisplayController::motion( const int& x, const int& y )
{
    // If a left button drag is in progress, translate the view
    if( m_left_drag )
    {
        const scalar& dx = x - m_last_x;
        const scalar& dy = y - m_last_y;
        m_last_x = x;
        m_last_y = y;
        translateView( dx, dy );
    }
    // If a right button drag is in progress, zoom the view
    if( m_right_drag )
    {
        const scalar& dx = x - m_last_x;
        const scalar& dy = y - m_last_y;
        m_last_x = x;
        m_last_y = y;
        zoomView( dx, dy );
    }
}

void TwoDimensionalDisplayController::translateView( const scalar& dx, const scalar& dy )
{
    const scalar& percent_x = dx/((scalar)m_window_width);
    const scalar& percent_y = dy/((scalar)m_window_height);
    const scalar& translate_x = percent_x*2.0*m_scale_factor*((scalar)m_window_width)/((scalar)m_window_height);
    const scalar& translate_y = percent_y*2.0*m_scale_factor;
    m_center_x -= translate_x;
    m_center_y += translate_y;
    reshape(m_window_width,m_window_height);
}

void TwoDimensionalDisplayController::zoomView( const scalar& dx, const scalar& dy )
{
    const scalar& percent_x = dx/((scalar)m_window_width);
    const scalar& percent_y = dy/((scalar)m_window_height);
    
    const scalar& scale = (std::fabs(percent_x) > std::fabs(percent_y)) ? -percent_x : percent_y;
    
    m_scale_factor += 2.0*scale;
    m_scale_factor = std::max(m_scale_factor,1.0e-4);
    
    reshape(m_window_width,m_window_height);
}

const int&  TwoDimensionalDisplayController::getWindowWidth() const
{
    return m_window_width;
}

const int&  TwoDimensionalDisplayController::getWindowHeight() const
{
    return m_window_height;
}

