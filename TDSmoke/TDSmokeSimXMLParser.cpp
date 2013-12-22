#include "TDSmokeSimXMLParser.h"

///////////////////////////////////////////////////////////////////////////////
// Main entry point

void TDSmokeSimXMLParser::loadExecutableSimulation( const std::string& file_name, SimulationEnsemble** sim_ensemble, scalar& max_time, int& velocity_pattern, std::string& description, std::string& scenetag, scalar &diff, scalar &visc, scalar &sigma, scalar &vf, scalar &vd, scalar &vg, bool &vg_enabled ) {
    // Load the xml document
    std::vector<char> xmlchars;
    rapidxml::xml_document<> doc;
    loadXMLFile( file_name, xmlchars, doc );
    
    // Attempt to locate the root node
    rapidxml::xml_node<>* node = doc.first_node("scene");
    if( node == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse xml scene file. Failed to locate root <scene> node. Exiting." << std::endl;
        exit(1);
    }
    
    // Determine what simulation type this is (particle, rigid body, etc)
    std::string simtype;
    loadSimulationType( node, simtype );
    
    // Parse common state
    loadMaxTime( node, max_time );
    loadSceneDescriptionString( node, description );
    loadSceneTag( node, scenetag );
    
    // Parse the user-requested simulation type. The default is a particle simulation.
    if( simtype == "td-smoke" ) {
        loadFluidRegionSize( node, sim_ensemble );
        loadVelocityPattern( node, velocity_pattern );
        loadDiffusion( node, diff );
        loadViscosity( node, visc );
        loadSmoothing( node, sigma );
        loadDrivingForceCoeff( node, vf );
        loadAttenuation( node, vd );
        loadGathering( node, vg, vg_enabled );
        loadMarker( node, sim_ensemble );
        loadTarget( node, sim_ensemble );
    }
    else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Invalid simtype '" << simtype << "' specified. Valid options are 'particle-system' and 'rigid-body'. Exiting." << std::endl;
        exit(1);
    }
}


///////////////////////////////////////////////////////////////////////////////
// Load XML file

void TDSmokeSimXMLParser::loadXMLFile( const std::string& filename, std::vector<char>& xmlchars, rapidxml::xml_document<>& doc ) {
    // Attempt to read the text from the user-specified xml file
    std::string filecontents;
    if( !loadTextFileIntoString(filename,filecontents) ) {
        std::cerr << "\033[31;1mERROR IN TWODSCENEXMLPARSER:\033[m XML scene file " << filename << ". Failed to read file." << std::endl;
        exit(1);
    }
    
    // Copy string into an array of characters for the xml parser
    for( int i = 0; i < (int) filecontents.size(); ++i ) xmlchars.push_back(filecontents[i]);
    xmlchars.push_back('\0');
    
    // Initialize the xml parser with the character vector
    doc.parse<0>(&xmlchars[0]);
}


///////////////////////////////////////////////////////////////////////////////
// Load text file into string

bool TDSmokeSimXMLParser::loadTextFileIntoString( const std::string& filename, std::string& filecontents ) {
    // Attempt to open the text file for reading
    std::ifstream textfile(filename.c_str(),std::ifstream::in);
    if(!textfile) return false;
    
    // Read the entire file into a single string
    std::string line;
    while(getline(textfile,line)) filecontents.append(line);
    
    textfile.close();
    
    return true;
}


///////////////////////////////////////////////////////////////////////////////
// Load scene tag

void TDSmokeSimXMLParser::loadSceneTag( rapidxml::xml_node<>* node, std::string& scenetag ) {
    assert( node != NULL );
    
    if( node->first_node("scenetag") ) {
        if( node->first_node("scenetag")->first_attribute("tag") ) {
            scenetag = node->first_node("scenetag")->first_attribute("tag")->value();
        }
        else {
            std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value of tag attribute for scenetag. Value must be string. Exiting." << std::endl;
            exit(1);
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
// Load max time

void TDSmokeSimXMLParser::loadMaxTime( rapidxml::xml_node<>* node, scalar& max_t ) {
    assert( node != NULL );
    
    // Attempt to locate the duraiton node
    rapidxml::xml_node<>* nd = node->first_node("duration");
    if( nd == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No duration specified. Exiting." << std::endl;
        exit(1);
    }
    
    // Attempt to load the duration value
    rapidxml::xml_attribute<>* timend = nd->first_attribute("time");
    if( timend == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No duration 'time' attribute specified. Exiting." << std::endl;
        exit(1);
    }
    
    max_t = std::numeric_limits<scalar>::signaling_NaN();
    if( !stringutils::extractFromString(std::string(timend->value()),max_t) ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'time' attribute for duration. Value must be numeric. Exiting." << std::endl;
        exit(1);
    }
}


///////////////////////////////////////////////////////////////////////////////
// Load scene description

void TDSmokeSimXMLParser::loadSceneDescriptionString( rapidxml::xml_node<>* node, std::string& description_string ) {
    assert( node != NULL );
    
    description_string = "No description specified.";
    
    // Attempt to locate the integrator node
    rapidxml::xml_node<>* nd = node->first_node("description");
    if( nd != NULL ) {
        rapidxml::xml_attribute<>* typend = nd->first_attribute("text");
        if( typend == NULL ) {
            std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No text attribute specified for description. Exiting." << std::endl;
            exit(1);
        }
        description_string = typend->value();
    }
}


///////////////////////////////////////////////////////////////////////////////
// Load fluid region size

void TDSmokeSimXMLParser::loadFluidRegionSize( rapidxml::xml_node<>* node, SimulationEnsemble** sim_ensemble ) {
    assert( node != NULL );
    
    int N;
    // Attempt to locate the fluid-region node
    rapidxml::xml_node<>* nd = node->first_node("fluid-region");
    if( nd == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No fluid-region specified. Exiting." << std::endl;
        exit(1);
    }
    
    // Attempt to load the fluid region size
    rapidxml::xml_attribute<>* Nnd = nd->first_attribute("n");
    if( Nnd == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No fluid-region 'n' attribute specified. Exiting." << std::endl;
        exit(1);
    }
    
    N = std::numeric_limits<int>::signaling_NaN();
    if( !stringutils::extractFromString(std::string(Nnd->value()),N) ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'n' attribute for fluid-region. Value must be numeric. Exiting." << std::endl;
        exit(1);
    }
    (*sim_ensemble) = new TDSmokeEnsemble(N);
}

///////////////////////////////////////////////////////////////////////////////
// Load velocity field pattern

void TDSmokeSimXMLParser::loadVelocityPattern( rapidxml::xml_node<>* node, int& velocity_pattern ) {
    assert( node != NULL );
    
    // Attempt to locate the duraiton node
    rapidxml::xml_node<>* nd = node->first_node("velocityfieldpattern");
    if( nd == NULL ) return;
    
    // Attempt to load the duration value
    rapidxml::xml_attribute<>* velocitynd = nd->first_attribute("code");
    if( velocitynd == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No pattern 'code' attribute specified. Exiting." << std::endl;
        exit(1);
    }
    
    velocity_pattern = std::numeric_limits<int>::signaling_NaN();
    if( !stringutils::extractFromString(std::string(velocitynd->value()), velocity_pattern) ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'code' attribute for velocity pattern. Value must be numeric. Exiting." << std::endl;
        exit(1);
    }
    
}


///////////////////////////////////////////////////////////////////////////////
// Load simulation type

void TDSmokeSimXMLParser::loadSimulationType( rapidxml::xml_node<>* node, std::string& simtype ) {
    assert( node != NULL );
    rapidxml::xml_node<>* nd = node->first_node("simtype");
    
    if( node->first_node("simtype") ) if( node->first_node("simtype")->first_attribute("type") ) simtype = node->first_node("simtype")->first_attribute("type")->value();
}


///////////////////////////////////////////////////////////////////////////////
// Load marker

void TDSmokeSimXMLParser::loadMarker( rapidxml::xml_node<>* node, SimulationEnsemble** sim_ensemble ) {
    assert( node != NULL );
    
    int marker = 0;
    for( rapidxml::xml_node<>* nd = node->first_node("marker"); nd; nd = nd->next_sibling("marker") ) {
        std::pair<int,int> markerpos;
        if( nd->first_attribute("i") ) {
            std::string attribute(nd->first_attribute("i")->value());
            if( !stringutils::extractFromString(attribute,markerpos.first) ) {
                std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value of i attribute for marker " << marker << ". Value must be integer. Exiting." << std::endl;
                exit(1);
            }
        }
        else {
            std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value of i attribute for marker " << marker << ". Exiting." << std::endl;
            exit(1);
        }
        
        if( nd->first_attribute("j") ) {
            std::string attribute(nd->first_attribute("j")->value());
            if( !stringutils::extractFromString(attribute,markerpos.second) ) {
                std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value of j attribute for marker " << marker << ". Value must be integer. Exiting." << std::endl;
                exit(1);
            }
        }
        else {
            std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value of j attribute for marker " << marker << ". Exiting." << std::endl;
            exit(1);
        }
        
        (*sim_ensemble)->setMarker( markerpos.first, markerpos.second );
        
        marker++;
    }
    (*sim_ensemble)->updateMarkerDensities();
}


///////////////////////////////////////////////////////////////////////////////
// Load target

void TDSmokeSimXMLParser::loadTarget( rapidxml::xml_node<>* node, SimulationEnsemble** sim_ensemble ) {
    assert( node != NULL );
    
    int target = 0;
    for( rapidxml::xml_node<>* nd = node->first_node("target"); nd; nd = nd->next_sibling("target") ) {
        std::pair<int,int> targetpos;
        if( nd->first_attribute("i") ) {
            std::string attribute(nd->first_attribute("i")->value());
            if( !stringutils::extractFromString(attribute,targetpos.first) ) {
                std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value of i attribute for target " << target << ". Value must be integer. Exiting." << std::endl;
                exit(1);
            }
        }
        else {
            std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value of i attribute for target " << target << ". Exiting." << std::endl;
            exit(1);
        }
        
        if( nd->first_attribute("j") ) {
            std::string attribute(nd->first_attribute("j")->value());
            if( !stringutils::extractFromString(attribute,targetpos.second) ) {
                std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value of j attribute for target " << target << ". Value must be integer. Exiting." << std::endl;
                exit(1);
            }
        }
        else {
            std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value of j attribute for target " << target << ". Exiting." << std::endl;
            exit(1);
        }
        
        (*sim_ensemble)->setTarget( targetpos.first, targetpos.second );
        
        target++;
    }
    
    (*sim_ensemble)->updateTargetDensities();
}


///////////////////////////////////////////////////////////////////////////////
// Load diffusion coefficient

void TDSmokeSimXMLParser::loadDiffusion( rapidxml::xml_node<>* node, scalar& diff ) {
    assert( node != NULL );
    
    // Attempt to locate the diffusion node
    rapidxml::xml_node<>* nd = node->first_node("diffusion");
    if( nd == NULL )
        return;
    
    // Attempt to load the diffusion value
    rapidxml::xml_attribute<>* diffnd = nd->first_attribute("diff");
    if( diffnd == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No diffusion 'diff' attribute specified. Exiting." << std::endl;
        exit(1);
    }
    
    diff = std::numeric_limits<scalar>::signaling_NaN();
    if( !stringutils::extractFromString(std::string(diffnd->value()),diff) ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'diff' attribute for diffusion. Value must be numeric. Exiting." << std::endl;
        exit(1);
    }
}


///////////////////////////////////////////////////////////////////////////////
// Load viscosity coefficient

void TDSmokeSimXMLParser::loadViscosity( rapidxml::xml_node<>* node, scalar& visc ) {
    assert( node != NULL );
    
    // Attempt to locate the viscosity node
    rapidxml::xml_node<>* nd = node->first_node("viscosity");
    if( nd == NULL ) return;
    
    // Attempt to load the viscosity value
    rapidxml::xml_attribute<>* viscnd = nd->first_attribute("visc");
    if( viscnd == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No viscosity 'visc' attribute specified. Exiting." << std::endl;
        exit(1);
    }
    
    visc = std::numeric_limits<scalar>::signaling_NaN();
    if( !stringutils::extractFromString(std::string(viscnd->value()),visc) ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'visc' attribute for viscosity. Value must be numeric. Exiting." << std::endl;
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Load smoothing parameter

void TDSmokeSimXMLParser::loadSmoothing( rapidxml::xml_node<>* node, scalar& sigma ) {
    assert( node != NULL );
    
    // Attempt to locate the smoothing node
    rapidxml::xml_node<>* nd = node->first_node("smoothing");
    if( nd == NULL ) {
        return;
    }
    
    // Attempt to load the smoothing value
    rapidxml::xml_attribute<>* sigmand = nd->first_attribute("sigma");
    if( sigmand == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No smoothing 'sigma' attribute specified. Exiting." << std::endl;
        exit(1);
    }
    
    sigma = std::numeric_limits<scalar>::signaling_NaN();
    if( !stringutils::extractFromString(std::string(sigmand->value()),sigma) ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'sigma' attribute for smoothing. Value must be numeric. Exiting." << std::endl;
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Load driving force coefficient

void TDSmokeSimXMLParser::loadDrivingForceCoeff( rapidxml::xml_node<>* node, scalar& vf ) {
    assert( node != NULL );
    
    // Attempt to locate the driving-force-coeff node
    rapidxml::xml_node<>* nd = node->first_node("drivingforce");
    if( nd == NULL ) return;
    
    // Attempt to load the driving force coefficient
    rapidxml::xml_attribute<>* vfnd = nd->first_attribute("vf");
    if( vfnd == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No drivingforce 'vf' attribute specified. Exiting." << std::endl;
        exit(1);
    }
    
    vf = std::numeric_limits<scalar>::signaling_NaN();
    if( !stringutils::extractFromString(std::string(vfnd->value()),vf) ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'vf' attribute for drivingforce. Value must be numeric. Exiting." << std::endl;
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Load attenuation coefficient

void TDSmokeSimXMLParser::loadAttenuation( rapidxml::xml_node<>* node, scalar& vd ) {
    assert( node != NULL );
    
    // Attempt to locate the attenuation node
    rapidxml::xml_node<>* nd = node->first_node("attenuation");
    if( nd == NULL ) return;
    
    // Attempt to load the attenuation value
    rapidxml::xml_attribute<>* vdnd = nd->first_attribute("vd");
    if( vdnd == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No attenuation 'vd' attribute specified. Exiting." << std::endl;
        exit(1);
    }
    
    vd = std::numeric_limits<scalar>::signaling_NaN();
    if( !stringutils::extractFromString(std::string(vdnd->value()),vd) ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'vd' attribute for attenuation. Value must be numeric. Exiting." << std::endl;
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Load gathering settings

void TDSmokeSimXMLParser::loadGathering( rapidxml::xml_node<>* node, scalar& vg, bool& vg_enabled ) {
    assert( node != NULL );
    
    // Attempt to locate the gathering node
    rapidxml::xml_node<>* nd = node->first_node("gathering");
    if( nd == NULL ) return;
    
    // Attempt to load the gathering coefficient
    rapidxml::xml_attribute<>* vgnd = nd->first_attribute("vg"); 
    if( vgnd == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No gathering 'vg' attribute specified. Exiting." << std::endl;
        exit(1);
    }
    
    vg = std::numeric_limits<scalar>::signaling_NaN();
    if( !stringutils::extractFromString(std::string(vgnd->value()),vg) ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'vg' attribute for gathering. Value must be numeric. Exiting." << std::endl;
        exit(1);
    }
    
    // Attempt to load the gathering enabled tag
    rapidxml::xml_attribute<>* enablednd = nd->first_attribute("enabled"); 
    if( enablednd == NULL ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No gathering 'enabled' attribute specified. Exiting." << std::endl;
        exit(1);
    }
    
    vg_enabled = std::numeric_limits<scalar>::signaling_NaN();
    if( !stringutils::extractFromString(std::string(enablednd->value()),vg_enabled) ) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'enabled' attribute for gathering. Value must be boolean. Exiting." << std::endl;
        exit(1);
    }
}

