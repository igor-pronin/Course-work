#ifndef ADGES_METHOD_H
#define ADGES_METHOD_H
#include <deque>
#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include "synonymous_types.h"
#include "base_types.h"
#include "Problem.h"
#include <ctime>
#include <cmath>
#include <string>
#include <algorithm>
using namespace std;

class TMethod {
    //размерность задачи
    int N;
    //количество функций ограничения
    int m;
    //глубина хранилища констант L
    int q;
    //тип выбора начальных точек
    char startCoverageType;
    //выбор глубины расширенной окрестности
    int r = 2;
    //вспомогательная переменная для создания расширенной окрестности
    int rcopy = r;
    //расширенная окрестность точки испытания
    set<TFuncsEvalutuonBlock* > T;
    //индексы начальных точек
    vector<vector<uchar>> indPoints;
    //нормированные координаты начальных точек
    vector<CoordinatesValues> points;
    //набор точек TFuncsEvalutuonBlock
    deque<TFuncsEvalutuonBlock> FuncsEvalutionBlocks;
    //набор рёбер TEdge
    deque<TEdge> EdgeBlocks;
    //набор точек, на основе которых формируются окрестности
    set<TFuncsEvalutuonBlock * > X;
    //набор точек, чьи окрестности нужно обновить после добавления новой точки
    set<TFuncsEvalutuonBlock * > X1;
    //класс решаемой задачи
    TProblem P;
    //дополнительный контейнер для формирования множества X при добавлении новой точки
    set<TFuncsEvalutuonBlock *> additionX;
    //отсортированный набор ребер
    deque<TEdge> GoodEdges;
    //количетсво поддерживаемых ребер
    int topEdges = 5;
    //нормированные координаты новой точки
    vector<double> coordinatesOfNewPoint;
    //индексы новой точки
    vector<uchar> indOfNewPoint;
    //максимальная константа Липшица
    double L = 4.0;
    //выбор характеристики ребра
    string edgeCharacteristic;
    //равномерное распределение
    bool uniformDistribution = false;

    
public:
    TMethod(TProblem & _P, uint lipshQueueDepth, char _startCoverageType, string _edgeCharacteristic) {
        N = _P.dimantion();
        m = _P.constraitsCount();
        q = lipshQueueDepth;
        startCoverageType = _startCoverageType;
        edgeCharacteristic = _edgeCharacteristic;
        P = _P;
        TFuncsEvalutuonBlock::setStaticFieldsValue(N, m, q);
        if (startCoverageType == 'b') {
            indPoints.resize(pow(2, N) + 1);
            points.resize(pow(2, N) + 1);
        }
        else {
            indPoints.resize(pow(3, N));
            points.resize(pow(3, N));
        }
        for (auto& g : indPoints) g.resize(N);
        for (auto& g : points) g.resize(N);
        startPoints(indPoints, points);
        FuncsEvalutionBlocks.resize(indPoints.size());
        for (uint i = 0; i < indPoints.size(); i++) {
            FuncsEvalutionBlocks[i] = (TFuncsEvalutuonBlock(points[i], indPoints[i], i));
            P.doEvaluation(FuncsEvalutionBlocks[i].f, FuncsEvalutionBlocks[i].grf, FuncsEvalutionBlocks[i].x);
            FuncsEvalutionBlocks[i].dx = P.coordinatesDecoding(FuncsEvalutionBlocks[i].x);
        }
        for (int i = 0; i < indPoints.size(); i++) {
            X.insert(& FuncsEvalutionBlocks[i]);
        }
        EdgeBlocks.resize(indPoints.size() * 2 * N);
        for (auto a : X) {
            neighborhood(*a);
        }
        indOfNewPoint.resize(2);
        createGoodEdges();
        if (edgeCharacteristic == "Piyavsky") {
            newPointByPiyavsky();
        }
        else if (edgeCharacteristic == "Length") {
            newPointByCenter();
        }
    }

    //сложение в двоичной системе счисления
    vector<uchar> addBin(vector<uchar>& A1, vector<uchar>& A2)
    {
        int i, i1, i2, k, p;
        static vector<uchar> R(N);
        for (int i = 0; i < N; i++)
            R[i] = 0;
        p = 0;
        i1 = N - 1;
        i2 = N - 1;
        for (i = N - 1; i >= 0; i--)
        {
            k = 0;
            if (i1 >= 0) k += A1[i1];
            if (i2 >= 0) k += A2[i2];
            k += p;
            R[i] = k % 2;
            p = k / 2;
            i1--;
            i2--;
        }
        return R;
    }

    //перевод числа из десятичной системы в троичную
    vector<uchar> ternarySystem(int i) {
        static vector<uchar> result(N);
        for (int i = 0; i < N; i++)
            result[i] = 0;
        int n = N - 1;
        do {
            result[n] = i % 3;
            i = i / 3;
            n--;
        } while (i >= 3);
        result[n] = i;
        return result;
    }


    //создание начальных точек
    void startPoints(vector<vector<uchar>> & indPoints, vector<CoordinatesValues>& points) {
        if (startCoverageType == 'b') {
           static vector<uchar> unit = { 0 , 1};
            for (int i = 1; i < indPoints.size() - 1; i++)
                indPoints[i] = addBin(indPoints[i - 1], unit);
            for (int i = 0; i < N; i++)
                indPoints[indPoints.size() - 1][i] = 2;
        }
        else {
            for (int i = 1; i < indPoints.size(); i++)
                indPoints[i] = ternarySystem(i);
        }
        for (int i = 0; i < indPoints.size(); i++)
            for (int j = 0; j < N; j++)
                if (indPoints[i][j] == 2)
                    points[i][j] = 0.5;
                else points[i][j] = indPoints[i][j];
    }

    //функция расстояния
    double distance(TFuncsEvalutuonBlock & xj, TFuncsEvalutuonBlock & xi) {
        double length = 0;
        for (int i = 0; i < N; i++)
            length += (xj.x[i] - xi.x[i]) * (xj.x[i] - xi.x[i]);
        return (sqrt(length) * alfa(xj.indX, xi.indX) * beta(xj, xi));
    }

    //функция альфа для функции расстояния
    double alfa(vector<uchar> & xj, vector<uchar> & xi) {
        double count = 0;
        for (int i = 0; i < N; i++)
            if (xj[i] == xi[i] && xj[i] != 2) count++;
        return (1.0 / (1.0 + count)); 
    }
    
    //функция бета для функции расстояния
    double beta(TFuncsEvalutuonBlock & xj, TFuncsEvalutuonBlock & xi) {
        if (xj.edges.size() == 0)
            return 1; 
        else {
            double max = 0;
            double option = 0;
            for (int i = 0; i < xj.edges.size(); i++) {
                option = cos(xj.x, FuncsEvalutionBlocks[EdgeBlocks[xj.edges[i]].numEndBlock].x, xi.x);
                if (option > max)
                    max = option;
            }
            //max = max > 1.0 ? 1.0 : max;
            return (1.0 / fabs(1.0 - max));
        }
    }

    double lengthBetweenPoints(TFuncsEvalutuonBlock & start, TFuncsEvalutuonBlock & finish) {
        double result = 0;
        for (int i = 0; i < N; i++) {
            double part = (finish.x[i] - start.x[i]) * (finish.x[i] - start.x[i]);
            result += part;
        }
        return sqrt(result);
    }

    //подсчет косинуса
    double cos(vector<TValue> & xj, vector<TValue> & ux, vector<TValue> & xi) {
        double result = 0;
        double result1 = 0;
        double result2 = 0;
        double result3 = 0;
        double result4 = 0;
        double result5 = 0;
     static vector<double> ij(N);
     static vector<double> uj(N);  
        for (int i = 0; i < N; i++) {
            ij[i] = xi[i] - xj[i];
            uj[i] = ux[i] - xj[i];
        }
        for (int i = 0; i < N; i++) {
            result += ij[i] * uj[i];
            result1 += ij[i] * ij[i];
            result2 += uj[i] * uj[i];
        }
        if (result < 0)
            return 0;
        result3 = sqrt(result1);     
        result4 = sqrt(result2); 
        result5 = result / (result3 * result4);
        return result5;
    }

    //создание окрестности точки
        void neighborhood(TFuncsEvalutuonBlock & point) {
            for (int i = 0; i < 2 * N; i++) {
                set <TFuncsEvalutuonBlock*> ::iterator it = X.begin();
                set <TFuncsEvalutuonBlock*> ::iterator it2 = X.end();
                set <TFuncsEvalutuonBlock*> ::iterator minIt;
                TValue ex = 100.0;
                TValue dis;
                TValue flag = 0;
                while (it != it2) {
                    for (auto a : point.edges)
                        if (FuncsEvalutionBlocks[EdgeBlocks[a].numEndBlock].x == (*(*it)).x)
                            flag = 1;
                    if (point.x == (*(*it)).x)
                        it++;
                    else  {
                        if (flag == 1)
                        it++;
                        else {
                            dis = distance(point, *(*it));
                            if (dis < ex) {
                                ex = dis;
                                minIt = it;
                            }
                            it++;
                            }
                    }
                    flag = 0;
                }
                TEdge edge = TEdge(point.numInFEB, (*(*minIt)).numInFEB);
                edge.length = lengthBetweenPoints(point, (*(*minIt)));
                if (edgeCharacteristic == "Piyavsky") {
                    edge.R = calculateR(point, (*(*minIt)));
                }
                else if (edgeCharacteristic == "Length") {
                    edge.R = edge.length * (-1.0);
                }
                int ind = point.numInFEB * 2 * N + i;
                EdgeBlocks[ind] = edge;
                point.edges.push_back(ind);
            }
        }

        void neighborhoodForNewPoint(TFuncsEvalutuonBlock& point) {
            for (int i = 0; i < 2 * N; i++) {
                set <TFuncsEvalutuonBlock*> ::iterator it = X1.begin();
                set <TFuncsEvalutuonBlock*> ::iterator it2 = X1.end();
                set <TFuncsEvalutuonBlock*> ::iterator minIt;
                bool x1IsEnded = false;
                TValue ex = 100.0;
                TValue dis;
                TValue flag = 0;
                while (it != it2) {
                    for (auto a : point.edges)
                        if (FuncsEvalutionBlocks[EdgeBlocks[a].numEndBlock].x == (*(*it)).x)
                            flag = 1;
                    if (point.x == (*(*it)).x)
                        it++;
                    else {
                        if (flag == 1)
                            it++;
                        else {
                            dis = distance(point, *(*it));
                            if (dis < ex) {
                                ex = dis;
                                minIt = it;
                            }
                            it++;
                        }
                    }
                    flag = 0;
                    if (it == it2 && x1IsEnded == false) {
                        it = additionX.begin();
                        it2 = additionX.end();
                        x1IsEnded = true;
                    }
                }
                TEdge edge = TEdge(point.numInFEB, (*(*minIt)).numInFEB);
                edge.length = lengthBetweenPoints(point, (*(*minIt)));
                if (edgeCharacteristic == "Piyavsky") {
                    edge.R = calculateR(point, (*(*minIt)));
                }
                else if (edgeCharacteristic == "Length") {
                    edge.R = edge.length * (-1.0);
                }
                int ind = point.numInFEB * 2 * N + i;
                EdgeBlocks[ind] = edge;
                point.edges.push_back(ind);
            }
        }

        //вспомогательная функция для построения набора точек, чьи окрестности подлежат
        //обновлению после добавления новой точки
        void extendedNeighborhood(TFuncsEvalutuonBlock & point) {
            for (int i = 0; i < point.edges.size(); i++) {
                T.insert(&FuncsEvalutionBlocks[EdgeBlocks[point.edges[i]].numEndBlock]);
                if (rcopy > 1) {
                    rcopy--;
                    extendedNeighborhood(FuncsEvalutionBlocks[EdgeBlocks[point.edges[i]].numEndBlock]);
                }
                if (i == point.edges.size() - 1) {
                    rcopy = r;
                }
            }
        }

        //создание новой точки и обновление определенных окрестностей
        void addPoint() {
            TFuncsEvalutuonBlock p = TFuncsEvalutuonBlock(coordinatesOfNewPoint, indOfNewPoint, FuncsEvalutionBlocks.size());
            P.doEvaluation(p.f, p.grf, p.x);
            p.dx = P.coordinatesDecoding(p.x);
            FuncsEvalutionBlocks.push_back(p);
            X1.insert(&FuncsEvalutionBlocks[FuncsEvalutionBlocks.size() - 1]);
            extendedNeighborhood(FuncsEvalutionBlocks[GoodEdges[0].numStartBlock]);
            X1.insert(&FuncsEvalutionBlocks[GoodEdges[0].numStartBlock]);
            for (auto a : T)
                X1.insert(a);
            T.clear();
            extendedNeighborhood(FuncsEvalutionBlocks[GoodEdges[0].numEndBlock]);
            X1.insert(&FuncsEvalutionBlocks[GoodEdges[0].numEndBlock]);
            for (auto a : T)
                X1.insert(a);
            T.clear();
            for (int i = 0; i < 2 * N; i++) {
                EdgeBlocks.push_back(TEdge());
            }
            for (auto a : X1) {
                for (auto b : a->edges)
                    additionX.insert(&FuncsEvalutionBlocks[EdgeBlocks[b].numEndBlock]);
                a->edges.clear();
                neighborhoodForNewPoint(*a);
                additionX.clear();
            }
            X1.clear();

            updateGoodEdges();
            if (edgeCharacteristic == "Piyavsky") {
                newPointByPiyavsky();
            }
            else if (edgeCharacteristic == "Length") {
                newPointByCenter();
            }
        }

        void createGoodEdges() {
            GoodEdges = EdgeBlocks;
            auto comp_func = [&](TEdge a, TEdge b) {return a.R < b.R; };
            sort(GoodEdges.begin(), GoodEdges.end(), comp_func);
            GoodEdges.erase(GoodEdges.begin() + topEdges, GoodEdges.end());
        }
        
        void deleteDisapperaredEdges() {
            bool edgeExist = false;
            deque<TEdge> ::iterator start = GoodEdges.begin();
            while (start != GoodEdges.end()) {
                for (auto b : FuncsEvalutionBlocks[start->numStartBlock].edges) {
                    if (EdgeBlocks[b].numEndBlock == start->numEndBlock) {
                        edgeExist = true;
                        EdgeBlocks[b].inGoodEdges = true;
                        break;
                    }
                }
                    if (!edgeExist) {
                        start = GoodEdges.erase(start);
                    }
                    else {
                        start++;
                    }
                    edgeExist = false;
            }
        }

        void checkNewEdges() {
            for (auto& a : EdgeBlocks) {
                if (!a.inGoodEdges) {
                    if (GoodEdges.size() < 5) {
                        GoodEdges.push_back(a);
                        auto comp_func = [&](TEdge a, TEdge b) {return a.R < b.R; };
                        sort(GoodEdges.begin(), GoodEdges.end(), comp_func);
                    }
                    else {
                        if (a.R < GoodEdges.back().R) {
                            deque<TEdge> ::iterator it = GoodEdges.end() - 2;
                            while (a.R < it->R && it != GoodEdges.begin()) {
                                it--;
                            }
                            if (a.R < GoodEdges.begin()->R) {
                                GoodEdges.insert(it, a);
                            }
                            else {
                                GoodEdges.insert(it + 1, a);
                            }
                            GoodEdges.pop_back();
                        }
                    }
                }
            }
        }

        void updateGoodEdges() {
            TEdge oldEdge = GoodEdges[0];
            deleteDisapperaredEdges();
            checkNewEdges();
            TEdge newEdge = GoodEdges[0];
            if (oldEdge.R == newEdge.R) {
                uniformDistribution = false;
            }
            else {
                uniformDistribution = true;
            }
        }

        //вывод в файл points.txt всех точек и их окрестностей после добавления новой точки
        void outputAllPoints() {
            ofstream fout;
            fout.open("TXT/points.txt");
            if (!fout.is_open())
            {
                cout << "Файл не может быть открыт или создан\n";
            }
            for (int i = 0; i < FuncsEvalutionBlocks.size(); i++) {
                for (int j = 0; j < N; j++)
                    fout << FuncsEvalutionBlocks[i].x[j] << " ";
                fout << endl;

                for (int k = 0; k < FuncsEvalutionBlocks[i].edges.size(); k++) {
                    for (int j = 0; j < N; j++)
                        fout << FuncsEvalutionBlocks[EdgeBlocks[FuncsEvalutionBlocks[i].edges[k]].numEndBlock].x[j] << " ";
                    fout << endl;
                }
            }
            fout.close();
        }

        void outputBestEdgeOnStep(int number) {
            string numberString = to_string(number);
            ofstream fout;
            fout.open("TXT/bestEdges/bestEdge" + numberString + ".txt");
            if (!fout.is_open())
            {
                cout << "Файл не может быть открыт или создан\n";
            }
            fout << GoodEdges[0].numStartBlock << endl;
            fout << GoodEdges[0].numEndBlock << endl;
            fout.close();
        }

        void outputEdgesOnStep(int number) {
            string numberString = to_string(number);
            ofstream fout;
            fout.open("TXT/edgesOnStep/edgeOnStep" + numberString + ".txt");
            if (!fout.is_open())
            {
                cout << "Файл не может быть открыт или создан\n";
            }
            for (int i = 0; i < FuncsEvalutionBlocks.size(); i++) {
                for (int j = 0; j < N; j++)
                    fout << FuncsEvalutionBlocks[i].x[j] << " ";
                fout << endl;

                for (int k = 0; k < FuncsEvalutionBlocks[i].edges.size(); k++) {
                    for (int j = 0; j < N; j++)
                        fout << FuncsEvalutionBlocks[EdgeBlocks[FuncsEvalutionBlocks[i].edges[k]].numEndBlock].x[j] << " ";
                    fout << endl;
                }
            }
            fout.close();
        }

        void outputBaseValues(int count, int N) {
            ofstream fout;
            fout.open("TXT/baseValues.txt");
            if (!fout.is_open())
            {
                cout << "Файл не может быть открыт или создан\n";
            }
            fout << count << endl;
            fout << N << endl;
            fout.close();
        }

        void outputNextEddedPoint() {
            ofstream fout;
            fout.open("TXT/nextAddedPoint.txt");
            if (!fout.is_open())
            {
                cout << "Файл не может быть открыт или создан\n";
            }
            fout << coordinatesOfNewPoint[0] << " ";
            fout << coordinatesOfNewPoint[1] << endl;
            fout << GoodEdges[0].numStartBlock << endl;
            fout << GoodEdges[0].numEndBlock << endl;
            for (auto edge : FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].edges) {
                fout << EdgeBlocks[edge].numEndBlock << endl;
            }
            for (auto edge : FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].edges) {
                fout << EdgeBlocks[edge].numEndBlock << endl;
            }

            fout.close();
        }

        void outputExtendedNeighborhoodForLastPoint() {
            ofstream fout;
            fout.open("TXT/extendedNeighborhoodForLastPoint.txt");
            if (!fout.is_open())
            {
                cout << "Файл не может быть открыт или создан\n";
            }
            extendedNeighborhood(FuncsEvalutionBlocks[GoodEdges[0].numStartBlock]);
            X1.insert(&FuncsEvalutionBlocks[GoodEdges[0].numStartBlock]);
            for (auto a : T)
                X1.insert(a);
            T.clear();
            extendedNeighborhood(FuncsEvalutionBlocks[GoodEdges[0].numEndBlock]);
            X1.insert(&FuncsEvalutionBlocks[GoodEdges[0].numEndBlock]);
            for (auto a : T)
                X1.insert(a);
            T.clear();
            fout << coordinatesOfNewPoint[0] << " ";
            fout << coordinatesOfNewPoint[1] << endl;
            for (auto point : X1) {
                fout << point->x[0] << " ";
                fout << point->x[1] << endl;
            }
            X1.clear();
        }


        //вычисление характеристики ребра методом Пиявского
        double calculateR(TFuncsEvalutuonBlock& start, TFuncsEvalutuonBlock& finish) {
            double dx = sqrt((start.dx[0] - finish.dx[0]) * (start.dx[0] - finish.dx[0]) +
                (start.dx[1] - finish.dx[1]) * (start.dx[1] - finish.dx[1]));
            return ((start.f[0] + finish.f[0]) / 2 - L *  dx / 2);
        }
        //кооридинаты новой точки методом Пиявского
        void newPointByPiyavsky() {
            coordinatesOfNewPoint.clear();
            double x, y;
            double rx, ry;
            x = FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[0] - FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[0];
            y = FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1] - FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1];

            if (x == 0) {
                rx = FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[0];
                ry = P.coordinateEncoding((FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].dx[1] + FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].dx[1]) / 2 -
                    (FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].f[0] - FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].f[0]) / (2 * L), 1);
                if (y > 0) {
                    if (ry > FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1] || ry < FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1]) {
                        cout << "вышли за границу Y1" << endl;
                        ry = (FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1] + FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1]) / 2;
                    }
                }
                else {
                    if (ry < FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1] || ry > FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1]) {
                        cout << "вышли за границу Y1" << endl;
                        ry = (FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1] + FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1]) / 2;
                    }
                }
            }
            else if (y == 0) {
                ry = FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1];
                rx = P.coordinateEncoding((FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].dx[0] + FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].dx[0]) / 2 -
                    (FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].f[0] - FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].f[0]) / (2 * L), 0);
                if (x > 0) {
                    if (rx > FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[0] || rx < FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[0]) {
                        cout << "вышли за границу X1" << endl;
                        rx = (FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[0] + FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[0]) / 2;
                    }
                }
                else {
                    if (rx < FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[0] || rx > FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[0]) {
                        cout << "вышли за границу X1" << endl;
                        rx = (FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[0] + FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[0]) / 2;
                    }
                }
            }
            else {
                ry = P.coordinateEncoding((FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].dx[1] + FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].dx[1]) / 2 -
                    (FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].f[0] - FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].f[0]) / (2 * L), 1);
                if (y > 0) {
                    if (ry > FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1] || ry < FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1]) {
                        cout << "вышли за границу Y2" << endl;
                        ry = (FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1] + FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1]) / 2;
                    }
                }
                    else {
                        if (ry < FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1] || ry > FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1]) {
                            cout << "вышли за границу Y2" << endl;
                            ry = (FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1] + FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1]) / 2;
                        }
                    }
                rx = ((ry - FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1]) * x / y) +
                    FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[0];
            }

            coordinatesOfNewPoint.push_back(rx);
            coordinatesOfNewPoint.push_back(ry);

            for (int j = 0; j < N; j++) {
                if (coordinatesOfNewPoint[j] == 0)
                    indOfNewPoint[j] = 0;
                else if (coordinatesOfNewPoint[j] == 1)
                    indOfNewPoint[j] = 1;
                else
                    indOfNewPoint[j] = 2;
            }
        }

        //координаты новой точки по середине отрезка
        void newPointByCenter() {
            coordinatesOfNewPoint.clear();
            double x, y;
            double rx, ry;
            x = FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[0] - FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[0];
            y = FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1] - FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1];

            if (x == 0) {
                rx = FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[0];
                ry = y / 2;
                ry += FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1];
            }
            else if (y == 0) {
                ry = FuncsEvalutionBlocks[GoodEdges[0].numStartBlock].x[1];
                rx = x / 2;
                rx += FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[0];
            }
            else {
                ry = y / 2;
                ry += FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[1];
                rx = x / 2;
                rx += FuncsEvalutionBlocks[GoodEdges[0].numEndBlock].x[0];
            }

            coordinatesOfNewPoint.push_back(rx);
            coordinatesOfNewPoint.push_back(ry);

            for (int j = 0; j < N; j++) {
                if (coordinatesOfNewPoint[j] == 0)
                    indOfNewPoint[j] = 0;
                else if (coordinatesOfNewPoint[j] == 1)
                    indOfNewPoint[j] = 1;
                else
                    indOfNewPoint[j] = 2;
            }
        }

        bool returnUniformDistribution() {
            return uniformDistribution;
        }

        double returnBestEdgeLength() {
            return GoodEdges[0].length;
        }

};

#endif