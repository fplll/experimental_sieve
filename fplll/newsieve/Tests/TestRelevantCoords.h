#ifndef TEST_RELEVANT_COORDS_H
#define TEST_RELEVANT_COORDS_H

#include "../DebugAll.h"
#include "../Typedefs.h"
#include <iostream>
#include <type_traits>
#include "../SieveUtility.h"
#include "../EMVApproximation.h"
#include "gmpxx.h"
#include <vector>

#include "../RelevantCoords.h"

bool test_relevant_coords()
{
  int dim = 20;
  using RelevantCoords = GaussSieve::RelevantCoordinates;
  
  RelevantCoords::get_instance(dim);
  
  int val = RelevantCoords::get_instance(dim).get_ij_value(2,1);
  std::cout << "[2,1]: " << RelevantCoords::get_instance(dim).get_ij_value(2,1) << std::endl;
  std::cout << "[2,1]: " << RelevantCoords::get_instance(dim).get_ij_value(2,1) << std::endl;
  std::cout << "[10,3]: " << RelevantCoords::get_instance(dim).get_ij_value(10,3) << std::endl;
  
  return true;
}

#endif