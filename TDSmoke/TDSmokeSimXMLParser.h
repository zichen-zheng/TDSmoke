#ifndef __TD_SMOKE_SIM_XML_PARSER_H__
#define __TD_SMOKE_SIM_XML_PARSER_H__

#include <Eigen/StdVector>

#include <iostream>
#include <fstream>
#include <limits>

#include "MathDefines.h"
#include "StringUtilities.h"
#include "SimulationEnsemble.h"
#include "TDSmokeEnsemble.h"

#include "rapidxml.hpp"

class TDSmokeSimXMLParser
{
public:
  
  void loadExecutableSimulation( const std::string& file_name, SimulationEnsemble** sim_ensemble, scalar& max_time, int& velocity_pattern, std::string& description, std::string& scenetag, scalar &diff, scalar &visc, scalar &sigma, scalar &vf, scalar &vd, scalar &vg, bool &vg_enabled );

private:
  void loadXMLFile( const std::string& filename, std::vector<char>& xmlchars, rapidxml::xml_document<>& doc );
  bool loadTextFileIntoString( const std::string& filename, std::string& filecontents );
  void loadSceneTag( rapidxml::xml_node<>* node, std::string& scenetag );
  void loadMaxTime( rapidxml::xml_node<>* node, scalar& max_t );
  void loadSceneDescriptionString( rapidxml::xml_node<>* node, std::string& description_string );
  void loadFluidRegionSize( rapidxml::xml_node<>* node, SimulationEnsemble** sim_ensemble );
  void loadVelocityPattern( rapidxml::xml_node<>* node, int& velocity_pattern );
  void loadSimulationType( rapidxml::xml_node<>* node, std::string& simtype );
  void loadViscosity( rapidxml::xml_node<>* node, scalar& visc );
  void loadDiffusion( rapidxml::xml_node<>* node, scalar& diff );
  void loadSmoothing( rapidxml::xml_node<>* node, scalar& sigma );
  void loadDrivingForceCoeff( rapidxml::xml_node<>* node, scalar& vf );
  void loadAttenuation( rapidxml::xml_node<>* node, scalar& vd );
  void loadGathering( rapidxml::xml_node<>* node, scalar& vg, bool& vg_enabled );
  void loadMarker( rapidxml::xml_node<>* node, SimulationEnsemble** sim_ensemble );
  void loadTarget( rapidxml::xml_node<>* node, SimulationEnsemble** sim_ensemble );

};

#endif
