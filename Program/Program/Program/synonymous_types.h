#ifndef TYPES_H
#define TYPES_H
#include <functional>
#include <vector>
using uint = unsigned int;
using uchar = unsigned char;
using TValue = double;
using CoordinateValue = double;
using FunctionValue = double;
using CoordinatesValues = std::vector<CoordinateValue>;
using FunctionsValues = std::vector<FunctionValue>;
using GradFunctionsValues = std::vector<FunctionsValues>;
using FunctionsCalculator =
std::function< int (FunctionsValues &, GradFunctionsValues &, const CoordinatesValues &)>;

#endif

