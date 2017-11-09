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
  const int dim = 20;
  
  using RelevantCoords = GaussSieve::RelevantCoordinates;
  
  GaussSieve::StaticInitializer<RelevantCoords> init1(dim);
  RelevantCoords matrix_of_rel_coo;
  
  
  std::cout << "matrix_of_rel_coo[3,2] = " << matrix_of_rel_coo.get_ij_value(3,2) << std::endl;
  std::cout << "matrix_of_rel_coo[1,2] = " << matrix_of_rel_coo.get_ij_value(1,2) << std::endl;
  
  matrix_of_rel_coo.print();
  
  /*
  RelevantCoords::get_instance(dim);
  
  int val = RelevantCoords::get_instance(dim).get_ij_value(2,1);
  std::cout << "[2,1]: " << RelevantCoords::get_instance(dim).get_ij_value(2,1) << std::endl;
  std::cout << "[2,1]: " << RelevantCoords::get_instance(dim).get_ij_value(2,1) << std::endl;
  std::cout << "[10,3]: " << RelevantCoords::get_instance(dim).get_ij_value(10,3) << std::endl;
  */
  return true;
}

#endif