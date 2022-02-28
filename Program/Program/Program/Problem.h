#ifndef PROBLEM_H
#define PROBLEM_H
#include "synonymous_types.h"

class TProblem {
    //размерность задачи
    uint Fdimanthion;
    //количество функций ограничения
    uint FconstraintsCount;
    //левая граница гиперинтервала
    CoordinatesValues Fa;
    //правая граница гиперинтервала
    CoordinatesValues Fb;
    //вектор декодированных координат
    CoordinatesValues decodedX;
    //хранение функции
    FunctionsCalculator Fp;

public:
    //возращение размерности
    uint dimantion() {return Fdimanthion;}
    //возращение количества ограничений
    uint constraitsCount() { return FconstraintsCount;}
    TProblem() {};
    TProblem(uint dimanthion, uint constraintsCount,
        const CoordinatesValues & a,const CoordinatesValues & b,
            FunctionsCalculator p) {
        Fdimanthion = dimanthion; FconstraintsCount = constraintsCount;
        Fa = a;  Fb = b; Fp = p;
        decodedX.resize(dimanthion);
    }
    //вспомогательная функция
    CoordinateValue coordinateDecoding(const CoordinateValue & xEncoded, int i)
    {
        return Fa[i] + xEncoded * (Fb[i] - Fa[i]);
    }
    //вспомогательная функция
    CoordinateValue coordinateEncoding(const CoordinateValue& xDecoded, int i)
    {
        return  (xDecoded - Fa[i])/ (Fb[i] - Fa[i]);
    }
    //пересчет нормированных (закодированных) координат в ненормированные
    CoordinatesValues& coordinatesDecoding(const CoordinatesValues& xEncoded)
    {
        for (size_t i(0); i < xEncoded.size(); ++i)
            decodedX[i] = coordinateDecoding(xEncoded[i], i);
        return decodedX; //зачем мы возвращаем ссылку?
    }
    //пересчет исходных координат в нормированные
    void coordinatesEncoding(const CoordinatesValues& xDecoded, 
                                      CoordinatesValues& xEncoded)
    {
        for (size_t i(0); i < xDecoded.size(); ++i)
            xEncoded[i] = coordinateEncoding(xDecoded[i], i);
    }
    //выполнение функции в декодированных координатах
    void doDecodedEvaluation(FunctionsValues& funcs, GradFunctionsValues& grads,
        const CoordinatesValues& xDecoded) const
    {  Fp(funcs,grads, xDecoded);
    }
    //выполнение функции в нормированных координатах
    void doEvaluation(FunctionsValues& funcs, GradFunctionsValues& grads,
        const CoordinatesValues& xEncoded)
    {   coordinatesDecoding(xEncoded);
        Fp(funcs, grads, decodedX);
    }
    //возврат левой границы гиперинтервала
  const  CoordinatesValues & leftBorders() const {return Fa;}
  //возврат правой границы гиперинтервала
  const  CoordinatesValues & rightBorders() const {return Fb;}
};

#endif