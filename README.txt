In a project root directory:
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

If you want to change cmake options (e.g. Debug/Release mode, export PNG, etc.), then enter the following command in build/ directory:
    $ ccmake ..

--------------------------------------------------------------

Special notes for OS X 10.9 Mavericks:

1. Update your Xcode to 5.0.1 for the latest version of SDK, which contains the necessary frameworks namely GLUT and OpenGL.

2. cmake -D OPENGL_INCLUDE_DIR=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers -D GLUT_INCLUDE_DIR=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/GLUT.framework/Versions/A/Headers ..

3. Fix the incompatibility in the FOSSSim code:
   a) Change the 176 line of (Your base folder)/include/eigen/Eigen/src/Core/Transpositions.h from "size" to "size()".
   b) Change the "bool operator<(const Interval &other)" in /FOSSSim/ContinuousTimeUtilities.h to "bool operator<(const Interval &other) const" (if applicable).


OR:

1. Use HomeBrew to bring gcc-4.2 back:
   $ brew install apple-gcc42

2. In /usr/bin/:
   $ sudo mv g++ g++-default
   $ sudo ln -s /usr/local/Cellar/apple-gcc42/4.2.1-xxxx.x/bin/g++-4.2 g++

3. Do similar stuffs to gcc, c++, and cpp

4. You don't need to change anything else to compile the source code!

