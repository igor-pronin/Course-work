#include "synonymous_types.h"
#include "base_types.h"
#include "Problem.h"
#include "adges_method.h"
#include "test_problem_functions.h"
#include "test_main.h"
using namespace std;

int main()
{
    int dimantion = 2, constraitsCount = 1;
    CoordinatesValues xLeft(dimantion), xReight(dimantion);
    xLeft = { -1.0, 0.0 };
    xReight = { 2.0, 3.0 };
    double epsilon = 0.02;
    TProblem P(dimantion, constraitsCount, xLeft, xReight, myFunctionMain);
    FunctionsValues funcs(P.constraitsCount() + 1u);
    GradFunctionsValues grads(P.constraitsCount() + 1u);
    for (auto& g : grads) g.resize(P.dimantion());
    CoordinatesValues xEncoded(P.dimantion());

    uint lipshQueueDepth = 4;
    char startCoverageType = 'q';
    //string edgeCharacteristic = "Length";
    string edgeCharacteristic = "Piyavsky";
    TMethod method(P, lipshQueueDepth, startCoverageType, edgeCharacteristic);
    int count = 300;

    for (int i = 0; i < count; i++) {
        if (method.returnBestEdgeLength() < epsilon) {
            count = i;
        }

        //method.outputBestEdgeOnStep(i);
        //method.outputEdgesOnStep(i);
        method.addPoint();

        if (i == count - 1) {
            if (!method.returnUniformDistribution()) {
                count++;
            }
            else {
                method.outputNextEddedPoint();
            }
        }
    }

    method.outputExtendedNeighborhoodForLastPoint();

    method.outputBaseValues(count, dimantion);
    method.outputBestEdgeOnStep(count);
    method.outputEdgesOnStep(count);

    return 0;
}
