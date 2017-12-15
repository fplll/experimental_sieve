#ifndef GAUSS_SIEVE_DEFAULT_INCLUDES_H
#define GAUSS_SIEVE_DEFAULT_INCLUDES_H

// We collect all the standard-library headers in this file.
// Since everything is one single big translation unit anyway, it does not make
// sense to do this fine-grained on a per-file level (and we would not notice if we
// forget something)
// We include this in every .h file.
// DebugAll.h is just for setting debug macros and we have to be sure it's included everywhere, so
// it goes here.
// Compat.h is just macros we need everywhere and standard library replacements, so it goes here.

#include "Compat.h"
#include "DebugAll.h"
#include "assert.h"
#include <array>
#include <atomic>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <exception>
#include <forward_list>
#include <fstream>
#include <iomanip>  // only to implement compute_statistics function; to delete
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <math.h>  //for log2, sqrt
#include <mutex>
#include <queue>
#include <random>
#include <stack>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#endif
