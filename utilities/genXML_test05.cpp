#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <time.h>
#include <stdlib.h>
using namespace std;

ofstream file;

const int N = 200;
const float sigma = 1.0;
const float vf = 2.0;
const float vd = 1.5;
const float vg = 0.0001;

//===============
// Main Function
//===============

int main() {
    file.open ("test05.xml");
    
    // XML file header (some definitions)
    file << "<scene>" << endl;
    file << "  <simtype type=\"td-smoke\"/>" << endl;
    file << "  <description text=\"Columbia University logo\"/>" << endl;
    file << "  <duration time=\"20.0\"/>" << endl;
    file << endl;
    file << "  <fluid-region n=\"" << N << "\"/>" << endl;
    file << "  <smoothing sigma=\"" << sigma << "\"/>" << endl;
    file << "  <drivingforce vf=\"" << vf << "\"/>" << endl;
    file << "  <attenuation vd=\"" << vd << "\"/>" << endl;
    file << "  <gathering vg=\"" << vg << "\" enabled=\"1\"/>" << endl;
    file << endl;
    
    ifstream markerFile, targetFile;
    markerFile.open("proc_img/cu_text_200.txt");
    targetFile.open("proc_img/cu_crown_200.txt");
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            bool m, t;
            markerFile >> m;
            targetFile >> t;
            if (m) file << "  <marker i=\"" << i << "\" j=\"" << j << "\"/>" << endl;
            if (t) file << "  <target i=\"" << i << "\" j=\"" << j << "\"/>" << endl;
        }
    }
    
    file << endl;
    
    // End of XML file
    file << "</scene>" << endl;
    
    markerFile.close();
    targetFile.close();
    file.close();
    return 0;
}

