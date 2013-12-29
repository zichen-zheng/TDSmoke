#include "MathDefines.h"
#include "MathUtilities.h"
#include <cmath>
#include <climits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include "TimingUtilities.h"
#include "StringUtilities.h"
#include "RenderingUtilities.h"
#include "YImage.h"
#include "TDSmokeEnsemble.h"
#include "TDSmokeSimXMLParser.h"


///////////////////////////////////////////////////////////////////////////////
// Contains the actual simulation, renderer, parser, and serializer

SimulationEnsemble* g_simulation_ensemble;


///////////////////////////////////////////////////////////////////////////////
// Rendering State

// If true, render the scene with OpenGL.
bool g_opengl_rendering_enabled = true;
// Used if 'capping' framerate for display purpouses. 0 means no cap.
scalar g_sec_per_frame = 0;
// Clock time last timestep was taken. Used when
scalar g_last_time = timingutils::seconds();


///////////////////////////////////////////////////////////////////////////////
// PNG Output State

std::string g_movie_dir = "pngs";
int g_steps_per_movie_frame = 0;


///////////////////////////////////////////////////////////////////////////////
// Parser state
std::string g_xml_scene_file;
std::string g_description;
std::string g_scene_tag = "";
bool g_xml_scene = false;
int g_velocity_pattern = 0;
scalar g_diff = DEFAULT_DIFF;
scalar g_visc = DEFAULT_VISC;
scalar g_sigma = DEFAULT_SIGMA;
scalar g_vf = DEFAULT_VF;
scalar g_vd = DEFAULT_VD;
scalar g_vg = DEFAULT_VG;
bool g_vg_enabled = DEFAULT_VG_ENABLED;


///////////////////////////////////////////////////////////////////////////////
// OpenGL window
int g_window;
int g_menu_id;
int g_submenu_id;
int g_value = 0;


///////////////////////////////////////////////////////////////////////////////
// Simulation state

// If run interactive mode, indicates if the simulation is paused or running.
bool g_paused = true;
// Timestep of the simulation
scalar g_dt = 0.01;
// Final (integer) step of the simulation. INT_MAX if no max.
int g_final_step = INT_MAX;
// Integer (current) step of simulation.
int g_current_step = 0;
int g_frame_count = 0;
int g_prev_time = 0;
int g_curr_time = 0;
scalar g_fps;
bool g_simulation_ran_to_completion = false;


#ifdef PNGOUT
// TODO: Preallocate space for for the image, use a faster method to read the buffer
void dumpPNG( const std::string& filename )
{
    YImage image;
    image.resize(g_simulation_ensemble->getWindowWidth(), g_simulation_ensemble->getWindowHeight());
    
    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glPixelStorei(GL_PACK_ROW_LENGTH, 0);
    glPixelStorei(GL_PACK_SKIP_ROWS, 0);
    glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
    glReadBuffer(GL_BACK);
    
    glFinish();
    glReadPixels(0, 0, g_simulation_ensemble->getWindowWidth(), g_simulation_ensemble->getWindowHeight(), GL_RGBA, GL_UNSIGNED_BYTE, image.data());
    image.flip();
    
    if( !image.save(filename.c_str()) )
    {
        std::cerr << outputmod::startred << "TDSmoke ERROR: " << outputmod::endred << "Failed to save png image " << filename << ". Exiting." << std::endl;
        exit(1);
    }
}
#endif

///////////////////////////////////////////////////////////////////////////////
// Simulation functions

//void miscOutputCallback();
//void sceneScriptingCallback();
//void dumpPNG(const std::string &filename);

void stepSystem()
{
    // If the user wants to generate a PNG movie
#ifdef PNGOUT
        std::stringstream oss;
        oss << g_movie_dir << "/frame" << std::setw(5) << std::setfill('0') << g_current_step << ".png";
        //std::cout << outputmod::startpink << "TDSmoke STATUS: " << outputmod::endpink << "Saving png to " << oss.str() << ". Exiting." << std::endl;
        dumpPNG(oss.str());
#endif
    
    // Determine if the simulation is complete.
    if( g_current_step >= g_final_step )
    {
        if( !g_opengl_rendering_enabled ) std::cout << std::endl;
        std::cout << outputmod::startpink << "TDSmoke STATUS: " << outputmod::endpink << "Simulation complete at time " << g_current_step*g_dt << ". Exiting." << std::endl;
        g_simulation_ran_to_completion = true;
        exit(0);
    }
    
    // Step the system forward in time
    g_simulation_ensemble->stepSystem(g_dt);
    g_current_step++;
}

void headlessSimLoop()
{
    if( g_final_step == INT_MAX )
    {
        std::cerr << outputmod::startred << "TDSmoke ERROR: " << outputmod::endred << " Must specify a final time in headless mode." << std::endl;
        exit(1);
    }
    
    scalar nextpercent = 0.02;
    std::cout << outputmod::startpink << "Progress: " << outputmod::endpink;
    for( int i = 0; i < 50; ++i ) std::cout << "-";
    std::cout << std::endl;
    std::cout << "          ";
    while( true )
    {
        scalar percent_done = ((scalar)g_current_step)/((scalar)g_final_step);
        if( percent_done >= nextpercent )
        {
            nextpercent += 0.02;
            std::cout << "." << std::flush;
        }
        stepSystem();
    }
}



///////////////////////////////////////////////////////////////////////////////
// Rendering and UI functions

void reshape( int w, int h )
{
    g_simulation_ensemble->reshape(w,h);
    assert( renderingutils::checkGLErrors() );
}

void renderBitmapString( const scalar& x, const scalar& y, const scalar& z, void *font, std::string s )
{
    glRasterPos3f(x, y, z);
    for( std::string::iterator i = s.begin(); i != s.end(); ++i )
    {
        char c = *i;
        glutBitmapCharacter(font, c);
    }
    
    assert( renderingutils::checkGLErrors() );
}

void setOrthographicProjection()
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    
    gluOrtho2D(0,g_simulation_ensemble->getWindowWidth(),0,g_simulation_ensemble->getWindowHeight());
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    assert( renderingutils::checkGLErrors() );
}

void menu(int num) {
    if (num == 0) {
        glutDestroyWindow(g_window);
        exit(0);
    }
    else {
        g_value = num;
    }
    glutPostRedisplay();
}

void createMenu() {
    g_submenu_id = glutCreateMenu(menu);
    glutAddMenuEntry("Black and White", 2);
    glutAddMenuEntry("Summer", 3);
    glutAddMenuEntry("Winter", 4);
    glutAddMenuEntry("Spring", 5);
    glutAddMenuEntry("Autumn", 6);
    g_menu_id = glutCreateMenu(menu);
    glutAddMenuEntry("Reset", 1);
    glutAddSubMenu("Color Map", g_submenu_id);
    glutAddMenuEntry("Quit", 0);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void computeFPS() {
    g_frame_count ++;
    g_curr_time = glutGet(GLUT_ELAPSED_TIME);
    int time_interval = g_curr_time - g_prev_time;
    if(time_interval > 1000) {
        // Calculate the number of frames per second
        g_fps = g_frame_count / (time_interval / 1000.0f);
        // Set time
        g_prev_time = g_curr_time;
        // Reset frame count
        g_frame_count = 0;
    }
}

void drawFPS() {
    setOrthographicProjection();
    glColor3d(1.0-g_simulation_ensemble->getBackgroundColor().r(),1.0-g_simulation_ensemble->getBackgroundColor().g(),1.0-g_simulation_ensemble->getBackgroundColor().b());
    
    std::string fpsStr;
    if (g_paused)
        fpsStr = "Paused";
    else
        fpsStr = stringutils::convertToString(g_fps,2) + " fps";
    renderBitmapString( 20, 20, 0.0, GLUT_BITMAP_HELVETICA_18, fpsStr );
    
    glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
    
    assert( renderingutils::checkGLErrors() );
}

void drawDescription() {
    setOrthographicProjection();
    glColor3d(1.0-g_simulation_ensemble->getBackgroundColor().r(),1.0-g_simulation_ensemble->getBackgroundColor().g(),1.0-g_simulation_ensemble->getBackgroundColor().b());
    
    std::string str(g_description);
    renderBitmapString( 20, 20, 0.0, GLUT_BITMAP_HELVETICA_18, str );
    
    glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
    
    assert( renderingutils::checkGLErrors() );
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    
    g_simulation_ensemble->display();
    
#ifdef PNGOUT
    drawDescription();
#else
    drawFPS();
#endif
    
    glutSwapBuffers();
    
    assert( renderingutils::checkGLErrors() );
}

//void centerCamera()
//{
//  renderingutils::Viewport view;
//
//  g_executable_simulation->computeCameraCenter(view);
//
//  scalar ratio = ((scalar)g_display_controller.getWindowHeight())/((scalar)g_display_controller.getWindowWidth());
//
//  view.size = 1.2*std::max(ratio*view.rx,view.ry);
//
//  g_display_controller.setCenterX(view.cx);
//  g_display_controller.setCenterY(view.cy);
//  g_display_controller.setScaleFactor(view.size);
//}

void keyboard( unsigned char key, int x, int y )
{
    // Handle general keyboard options shared across simulations
    
    // Exit on escape key or q
    if( key == 27 || key == 'q' || key == 'Q' )
    {
        exit(0);
    }
    // Step the system forward one timestep
    else if( key == 's' || key =='S' )
    {
        stepSystem();
        glutPostRedisplay();
    }
    // Space bar pauses and unpauses
    else if( key == ' ' )
    {
        g_paused = !g_paused;
    }
    //  else if( key == 'c' || key == 'C' )
    //  {
    //    centerCamera();
    //    g_display_controller.reshape(g_display_controller.getWindowWidth(),g_display_controller.getWindowHeight());
    //    glutPostRedisplay();
    //  }
#ifdef PNGOUT
    // Save a png screenshot
    else if( key == 'i' || key == 'I' )
    {
        std::cout << outputmod::startpink << "TDSmoke STATUS: " << outputmod::endpink << "Saving screenshot as 'output.png'." << std::endl;
        dumpPNG("output.png");
    }
#endif
    
    // Allow the simulation to respond to keyboard input on its own
    g_simulation_ensemble->keyboard(key,x,y);
    
    glutPostRedisplay();
    assert( renderingutils::checkGLErrors() );
}

// Proccess 'special' keys
void special( int key, int x, int y )
{
    g_simulation_ensemble->special( key, x, y );
    glutPostRedisplay();
    assert( renderingutils::checkGLErrors() );
}

void mouse( int button, int state, int x, int y )
{
    if ( !g_xml_scene || button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON )
    {
        g_simulation_ensemble->mouse( button, state, x, y );
        glutPostRedisplay();
        assert( renderingutils::checkGLErrors() );
    }
}

void motion( int x, int y )
{
    g_simulation_ensemble->motion( x, y );
    glutPostRedisplay();
    assert( renderingutils::checkGLErrors() );
}

void idle()
{
    // Trigger the next timestep if the time since the last step exceeds the threshold
    const scalar& current_time = timingutils::seconds();
    if( !g_paused && current_time-g_last_time >= g_sec_per_frame )
    {
        computeFPS();
        g_last_time = current_time;
        stepSystem();
        glutPostRedisplay();
    }
    
    assert( renderingutils::checkGLErrors() );
}

void initializeOpenGLandGLUT( int argc, char** argv )
{
    // Initialize GLUT
    glutInit(&argc,argv);
    glutInitDisplayMode(g_simulation_ensemble->getGlutDisplayMode());
    glutInitWindowSize(g_simulation_ensemble->getWindowWidth(),g_simulation_ensemble->getWindowHeight());
    
    // Position the window in the center of the screen
    int scnwdh = glutGet(GLUT_SCREEN_WIDTH);
    int scnhgt = glutGet(GLUT_SCREEN_HEIGHT);
    if( scnwdh != 0 && scnhgt != 0 )
    {
        assert( g_simulation_ensemble->getWindowWidth() < scnwdh );
        int strtx = (scnwdh-g_simulation_ensemble->getWindowWidth())/2;
        assert( g_simulation_ensemble->getWindowHeight() < scnhgt );
        int strty = (scnhgt-g_simulation_ensemble->getWindowHeight())/2;
        glutInitWindowPosition(strtx,strty);
    }
    
    g_window = glutCreateWindow("Target-Driven Smoke Simulation");
    createMenu();
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutIdleFunc(idle);
    
    // Initialize OpenGL
    g_simulation_ensemble->initializeOpenGL();
    reshape(g_simulation_ensemble->getWindowWidth(),g_simulation_ensemble->getWindowHeight());
    
    assert( renderingutils::checkGLErrors() );
}


///////////////////////////////////////////////////////////////////////////////
// Parser functions

void loadScene( const std::string& file_name )
{
    // Maximum time in the simulation to run for. This has nothing to do with run time, cpu time, etc. This is time in the 'virtual world'.
    scalar max_time;
    // Maximum frequency, in wall clock time, to execute the simulation for. This serves as a cap for simulations that run too fast to see a solution.
    scalar steps_per_sec_cap = 60.0;
    // Contains the center and 'scale factor' of the view
    //renderingutils::Viewport view;
    
    // Load the simulation and pieces of rendring and UI state
    TDSmokeSimXMLParser xml_scene_parser;
    xml_scene_parser.loadExecutableSimulation( file_name, &g_simulation_ensemble,
                                              max_time, g_velocity_pattern, g_description, g_scene_tag, g_diff, g_visc, g_sigma, g_vf, g_vd, g_vg, g_vg_enabled );
    
    
    g_fps = steps_per_sec_cap;
    // To cap the framerate, compute the minimum time a single timestep should take
    g_sec_per_frame = 1.0/steps_per_sec_cap;
    // Integer number of timesteps to take
    g_final_step = ceil(max_time/g_dt);
    // We begin at the 0th timestep
    //g_current_step = 0;
}

void parseCommandLine( int argc, char** argv )
{
    try
    {
        TCLAP::CmdLine cmd("TDSmoke");
        
        // Directory to save png movie frames too
        TCLAP::ValueArg<std::string> movie("m", "movie", "Directory to save png movie too", false, "", "string", cmd);
        
        // Frequency (in simulation time) at which to save frames
        TCLAP::ValueArg<scalar> moviefreq("f", "moviefrequency", "Frequency (in simulation time) at which to save frames", false, -1.0, "scalar", cmd);
        
        // Begin the scene paused or running
        TCLAP::ValueArg<bool> paused("p", "paused", "Begin the simulation paused if 1, running if 0", false, true, "boolean", cmd);
        
        // Run the simulation with rendering enabled or disabled
        TCLAP::ValueArg<bool> display("d", "display", "Run the simulation with display enabled if 1, without if 0", false, true, "boolean", cmd);
        
        // Xml scene file
        TCLAP::ValueArg<std::string> scene("s", "scene", "Xml scene file", false, "", "string", cmd);
        
        //    // Save svgs to a movie directory
        //    TCLAP::ValueArg<std::string> movie("m", "moviedir", "Directory to output svg screenshot to", false, "", "string", cmd);
        
        cmd.parse(argc, argv);
        
        if (scene.isSet())
        {
            g_xml_scene_file = scene.getValue();
            g_xml_scene = true;
        }
        
        g_paused = paused.getValue();
        g_opengl_rendering_enabled = display.getValue();
        
        if( movie.isSet() )
        {
            g_movie_dir = movie.getValue();
            scalar fps = 1.0/g_dt;
            if( moviefreq.isSet() )
            {
                fps = moviefreq.getValue();
            }
            assert( fps > 0.0 );
            assert( fps*g_dt <= 1.0 );
            
            if( !mathutils::approxEqual(fmod(1.0/g_dt,fps),0.0,1.0e-9) )
            {
                std::cerr << outputmod::startred << "TDSmoke WARNING: " << outputmod::endred << "Time between frames and timestep not integer multiple." << std::endl;
            }
            
            g_steps_per_movie_frame = mathutils::round(1.0/(fps*g_dt));
        }
    }
    catch( TCLAP::ArgException& e )
    {
        std::cerr << outputmod::startred << "TDSmoke ERROR: " << outputmod::endred << " Parse error from TCLAP: " << e.what() << std::endl;
        exit(1);
    }
}



///////////////////////////////////////////////////////////////////////////////
// Various support functions

void cleanupAtExit()
{
    if( g_simulation_ensemble != NULL )
    {
        delete g_simulation_ensemble;
        g_simulation_ensemble = NULL;
    }
}

int main( int argc, char** argv )
{
    srand(time(0));

    #ifdef PNGOUT
    system("rm pngs/*");
    #endif

    // Parse command line arguments
    parseCommandLine( argc, argv );
    
    // Function to cleanup at progarm exit
    atexit(cleanupAtExit);
    
    // If requestedm load xml scene parameters
    if (g_xml_scene)
    {
        //g_simulation_ensemble = new TDSmokeEnsemble(100);
        loadScene( g_xml_scene_file );
        
        g_simulation_ensemble->setVelocityPattern( g_velocity_pattern );
        g_simulation_ensemble->setDiffusion( g_diff );
        g_simulation_ensemble->setViscosity( g_visc );
        g_simulation_ensemble->setSmoothing( g_sigma );
        g_simulation_ensemble->setDrivingForceCoeff( g_vf );
        g_simulation_ensemble->setAttenuation( g_vd );
        g_simulation_ensemble->setGathering( g_vg );
        g_simulation_ensemble->setGatheringEnabled( g_vg_enabled );
    }
    else
    {
        g_simulation_ensemble = new TDSmokeEnsemble(100);
    }
    
    // Initialization for OpenGL and GLUT
    if( g_opengl_rendering_enabled ) initializeOpenGLandGLUT(argc,argv);
    
    std::cout << outputmod::startgreen;
    std::cout << "+----------------------------------+" << std::endl;
    std::cout << "|  Target-Driven Smoke Simulation  |" << std::endl;
    std::cout << "+----------------------------------+" << std::endl;
    std::cout << outputmod::endgreen;
    
    // TODO: Print build info here
#ifdef CMAKE_BUILD_TYPE
    std::cout << outputmod::startblue << "Build type: " << outputmod::endblue << CMAKE_BUILD_TYPE << std::endl;
#endif
    std::cout << outputmod::startblue << "Rendering: " << outputmod::endblue << (g_opengl_rendering_enabled?"Enabled":"Disabled") << std::endl;
    if( g_opengl_rendering_enabled && g_movie_dir != "" )
    {
        std::cout << outputmod::startblue << "Movie output directory: " << outputmod::endblue << g_movie_dir << std::endl;
        std::cout << outputmod::startblue << "Movie output FPS: " << outputmod::endblue << (1.0/(g_dt*((scalar)g_steps_per_movie_frame))) << std::endl;
    }
    
    if( g_xml_scene )
    {
        std::cout << outputmod::startblue << "Scene: " << outputmod::endblue << g_xml_scene_file << std::endl;
        std::cout << outputmod::startblue << "Description: " << outputmod::endblue << g_description << std::endl;
    }
    else
    {
        g_simulation_ensemble->initMarkerAndTarget();
    }
    
    std::cout << outputmod::startblue << "Smoothing parameter: " << outputmod::endblue << g_sigma << std::endl;
    std::cout << outputmod::startblue << "Driving force coefficient: " << outputmod::endblue << g_vf << std::endl;
    std::cout << outputmod::startblue << "Attenuation coefficient: " << outputmod::endblue << g_vd << std::endl;
    std::cout << outputmod::startblue << "Gathering coefficient: " << outputmod::endblue << g_vg << std::endl;
    std::cout << outputmod::startblue << "Gathering: " << outputmod::endblue;
    if (g_vg_enabled) std::cout << "Enabled" << std::endl;
    else std::cout << "Disabled" << std::endl;
    
    if( g_opengl_rendering_enabled ) glutMainLoop();
    else headlessSimLoop();
    
    return 0;
}
