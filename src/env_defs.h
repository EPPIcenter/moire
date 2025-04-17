#pragma once

#ifndef ENV_DEFS_H_
#define ENV_DEFS_H_

#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
#include <Rinterface.h>
#endif

#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
#include <execution>
#include "RcppParallel.h"
#define HAS_EXECUTION 1
#endif

#endif