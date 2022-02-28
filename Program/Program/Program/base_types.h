#ifndef BASETYPES_H
#define BASETYPES_H
#include <deque>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include "synonymous_types.h"
using namespace std;

struct TFuncsEvalutuonBlock
{
    //размерность задачи
   static uint N;
   //количество функций ограничения
   static uint m;
   //глубина хранилища констант L
   static uint q;
   //координаты точки
    vector<CoordinateValue> x;
    //декодированные координаты точки
    vector<double> dx;
    //индексы координат точки из множества {0;1;2}
    vector<uchar> indX;
    //значения основной функции и функций ограничений
    vector<FunctionValue> f;
    //значения градиентов
    vector<vector<FunctionValue>> grf;
    //номера ребер связей в блоке EdgeBlocks класса Method
    vector<size_t> edges;
    //хранилище констант L
    vector<vector<TValue>> locL;
    //номер точки в блоке FuncsEvalutionBlocks класса Method
    int numInFEB;
    static void setStaticFieldsValue(int _N, int _m, int _q)
    {N = _N; m = _m; q = _q;}
    TFuncsEvalutuonBlock() {};
    TFuncsEvalutuonBlock(vector<CoordinateValue> & point, 
        vector<uchar> & indPoint, uint _numInFEB) {

        x.resize(N);                          
        indX.resize(N);
        f.resize(m + 1);
        grf.resize(m + 1);
        for (auto& g : grf) g.resize(N);
        locL.resize(m + 1);
        for (auto& g : locL) g.resize(q);
        for (int i = 0; i < N; i++) {
            x[i] = point[i];
            indX[i] = indPoint[i];
        }
     numInFEB = _numInFEB; 
    }
};

struct TEdge
{
    //номер точки в блоке FuncsEvalutionBlocks класса Method, из которой исходит ребро
    int numStartBlock;
    //номер точки в блоке FuncsEvalutionBlocks класса Method, которым оканчивается ребро
    int numEndBlock;
    ////номер ребра в массиве edges породившего блок
    //int numEdgeInStartBlock;
    //номер ребра в блоке GoodEdges класса Method
    int numInGoodEdge;
    //длина ребра
    TValue length;
    //приоритет ребра
    TValue R;
    //смещение точки-кандидата на измерение от начала ребра
    TValue S;
    //
    bool inGoodEdges = false;
    TEdge() {};
    TEdge(int _numStartBlock, int _numEndBlock) {
        numStartBlock = _numStartBlock;
        numEndBlock = _numEndBlock;
        //numEdgeInStartBlock = 0;
        numInGoodEdge = 0;
        length = 0;
        R = 0;
        S = 0;
    }
    bool operator == (TEdge e) {
        if ((numStartBlock == e.numStartBlock && numEndBlock == e.numEndBlock) || (numStartBlock == e.numEndBlock && numEndBlock == e.numStartBlock))
            return true;
        else return false;
    }
};

#endif