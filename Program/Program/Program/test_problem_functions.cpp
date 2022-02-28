#include "synonymous_types.h"
#include "test_problem_functions.h"
#include <algorithm>

double myFunction1(const CoordinatesValues& x) {
    return 1 * sqrt((x[0] - 1.5) * (x[0] - 1.5) + (x[1] - 0.5) * (x[1] - 0.5)) + 0;
}

double myFunction2(const CoordinatesValues& x) {
    return 0.5 * sqrt((x[0] - 1.5) * (x[0] - 1.5) + (x[1] - 2.5) * (x[1] - 2.5)) + 0.5;
}

double myFunction3(const CoordinatesValues& x) {
    return 0.5 * sqrt((x[0] + 0.5) * (x[0] + 0.5) + (x[1] - 2.5) * (x[1] - 2.5)) + 0.25;
}

int myFunctionMain(FunctionsValues& funcs, GradFunctionsValues& grads,
    const CoordinatesValues& x) {
    funcs[0] = std::min({ myFunction1(x), myFunction2(x), myFunction3(x) });
    return 0;
}



//{ for (int j = 0; j < funcs.size(); j++)
//    switch (j)
//    {
//    case 0:
//        funcs[j] = x[0] * x[0] + 5 * x[1] * x[1];
//        grads[j][0] = 2 * x[0];
//        grads[j][1] = 10 * x[1];
//        break;
//    case 1:
//        funcs[j] = sin(x[0]) * x[1] + x[0] * x[1] * x[1];
//        grads[j][0] = cos(x[0]) * x[1] + x[1] * x[1];
//        grads[j][1] = sin(x[0]) + 2 * x[0] * x[1];
//        break;
//    default:
//        for (auto& gr : grads) for (auto& d : gr) d = 0;
//        for (auto& f : funcs) f = 0;
//    }
//    return 0;
//}
