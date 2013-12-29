======================================================================

                   Target-Driven Smoke Simulation

                      written by Zichen Zheng

======================================================================


---------
  ABOUT
---------

This is a 2D implementation of a SIGGRAPH 2004 paper "Target-Driven Smoke Animation" by et al.


------------------------
  WHAT'S IN THE SOURCE
------------------------

The following lists what is included in each directory:
TDSmoke:   source code
assets:    test files (in XML)
cmake:     CMake configuration files
include:   source code of required libraries (Eigen, RapidXML, TCLAP)
utilities: help programs (especially for creating test files)


------------------
  HOW TO COMPILE
------------------

In a project root directory:
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

If you want to change cmake options (e.g. Debug/Release mode, export PNG, etc.), then enter the following command in build/ directory:
    $ ccmake ..


--------------------------
  FOR OS X 10.9 MAVERICS
--------------------------

1. Update your Xcode to 5.0.1 for the latest version of SDK, which contains the necessary frameworks namely GLUT and OpenGL.

2. $ cmake -D OPENGL_INCLUDE_DIR=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers -D GLUT_INCLUDE_DIR=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/GLUT.framework/Versions/A/Headers ..

3. Fix the incompatibility in the TDSmoke code:
   Change the 176 line of (Your base folder)/include/eigen/Eigen/src/Core/Transpositions.h from "size" to "size()".


OR:

1. Use HomeBrew to bring gcc-4.2 back:
   $ brew install apple-gcc42

2. In /usr/bin/:
   $ sudo mv g++ g++-default
   $ sudo ln -s /usr/local/Cellar/apple-gcc42/4.2.1-xxxx.x/bin/g++-4.2 g++

3. Do similar stuffs to gcc, c++, and cpp

4. You don't need to change anything else to compile the source code!


--------------
  REFERENCES
--------------

Raanan Fattal and Dani Lischinski. 2004. Target-driven smoke animation. ACM Trans. Graph. 23, 3 (August 2004), 441-448.

