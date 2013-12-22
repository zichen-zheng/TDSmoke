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

const int N = 80;
const float sigma = 1.0;
const float vf = 2.0;
const float vd = 1.5;
const float vg = 0.0001;

//===============
// Main Function
//===============

int main() {
    file.open ("test01.xml");
    
    // XML file header (some definitions)
    file << "<scene>" << endl;
    file << "  <simtype type=\"td-smoke\"/>" << endl;
    file << "  <description text=\"Target-driven smoke: 'X' becomes 'O'\"/>" << endl;
    file << "  <duration time=\"20.0\"/>" << endl;
    file << endl;
    file << "  <fluid-region n=\"" << N << "\"/>" << endl;
    file << "  <smoothing sigma=\"" << sigma << "\"/>" << endl;
    file << "  <drivingforce vf=\"" << vf << "\"/>" << endl;
    file << "  <attenuation vd=\"" << vd << "\"/>" << endl;
    file << "  <gathering vg=\"" << vg << "\" enabled=\"1\"/>" << endl;
    file << endl;
    
    // marker
    for (int i = N/20; i < N-N/20; i++) {
        file << "  <marker i=\"" << i << "\" j=\"" << i << "\"/>" << endl;
        file << "  <marker i=\"" << i << "\" j=\"" << N-i << "\"/>" << endl;
    }
    
    // target
    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            float dist_to_center = pow(i-N/2,2) + pow(j-N/2,2);
            dist_to_center = sqrt(dist_to_center);
            if (fabs(dist_to_center - (float)N/4) < 1) {
                file << "  <target i=\"" << i << "\" j=\"" << j << "\"/>" << endl;
            }
        }
    }
    
    file << endl;
    
    // End of XML file
    file << "</scene>" << endl;
    
    file.close();
    return 0;
}

