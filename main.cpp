#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>   // rand
#include <math.h>    // sqrt, pow
#include <limits>    // numeric_limits<float>::max()
#include <ctime>
#include <algorithm>
using namespace std;

class City
{
public:
    City(int cityNum, float x, float y)
    {
        this->cityNum = cityNum;
        this->x = x;
        this->y = y;
    }

    City() { }

    int cityNum;
    float x;
    float y;
};

class TabuSearch
{
private:
    int tabuListLength;
    int numCities;
    int numIterations;
    int numMaxTries;
    int *globalSolution;
    City **cities;
    int **tabuList;
    bool initialized;
public:
    TabuSearch(int numIterations, int numMaxTries); // initialization
    ~TabuSearch();
    void initFromFile(string fileName);
    void initFromStream();
    int* genGreedySolution(int *solution, bool initialize); // generates greedy entry solution
    int* genRandomSolution(int *solution, bool initialize);
    float calcAllDist(int *cityLabel); // returns distance of the whole route
    float calcEdgeDist(City *&a, City *&b) const; // returns distance between two cities
    float calcNewDist(int *solution, int firstPair, int secondPair, float oldDist);
    void calcTabuIters(int *localSolution, int i, int &tabuIFirst, int &tabuISec);
    void calcTabuIters2(int *localSolution, int i, int j, int &tabuIFirst, int &tabuISec);
    float checkAllPermutations();
    void printBestSolution();
    void tabuSearch();
};

void TabuSearch::printBestSolution()
{
    for(int i=0; i<numCities+1; ++i)
        cout << globalSolution[i] << " ";
    cout << endl;
}

TabuSearch::TabuSearch(int numIterations, int numMaxTries)
{
    this->numIterations = numIterations;
    this->numMaxTries = numMaxTries;
    initialized = false;
}

TabuSearch::~TabuSearch()
{
    if(initialized)
    {
        delete []globalSolution;
        for(int i=0; i<numCities-1; ++i)
            delete tabuList[i];
        delete []tabuList;
        delete []cities;
    }
}

void TabuSearch::initFromFile(string fileName)
{
    ifstream file(fileName.c_str());
    if(file)
    {
        file >> numCities;
        tabuListLength = 3*numCities;
        if(tabuListLength <= 7) tabuListLength = 3;
        globalSolution = new int[numCities+1];
        cities = new City*[numCities+1];
        tabuList = new int *[numCities-1];
        for(int i=0; i<numCities-1; ++i)
        {
            tabuList[i] = new int[numCities-1-i];
            for(int j=0; j<numCities-1-i; ++j)
                tabuList[i][j] = 0;
        }

        int cityNum;
        float x, y;
        for(int i=0; i<numCities; ++i)
        {
           file >> cityNum >> x >> y;
           cities[i] = new City(cityNum, x, y);
        }
        cities[numCities] = cities[0];
        initialized = true;

        /*cout << "num cities: " << cityNum << endl;
        for(int i=0; i<numCities+1; ++i)
        {
            cout << cities[i]->cityNum << " " << cities[i]->x << " " << cities[i]->y << endl;
        } */
    }
    else
        cerr << "cant open file" << endl;
}

void TabuSearch::initFromStream()
{
    cin >> numCities;
    tabuListLength = 3*numCities;
    if(tabuListLength <= 7) tabuListLength = 3;
    globalSolution = new int[numCities+1];
    cities = new City*[numCities+1];
    tabuList = new int *[numCities-1];
    for(int i=0; i<numCities-1; ++i)
    {
        tabuList[i] = new int[numCities-1-i];
        for(int j=0; j<numCities-1-i; ++j)
            tabuList[i][j] = 0;
    }

    int cityNum;
    float x, y;
    for(int i=0; i<numCities; ++i)
    {
       cin >> cityNum >> x >> y;
       cities[i] = new City(cityNum, x, y);
    }
    cities[numCities] = cities[0];
}

float TabuSearch::calcEdgeDist(City *&a, City *&b) const
{
    return sqrt(pow((a->x - b->x), 2) + pow((a->y - b->y), 2));
}

float TabuSearch::calcAllDist(int *cityLabel)
{
    float dist = 0.0;
    for(int i=0; i<numCities; ++i) // go through all edges, city's length is numCities+1
    {
        dist += calcEdgeDist(cities[ cityLabel[i] ], cities[ cityLabel[i+1] ]);
    }
    return dist;
}

// generate greedy solution based on Cities **cities
int* TabuSearch::genGreedySolution(int *solution, bool initialize)
{
    if(initialize)
    {
        for(int i=0; i<numCities; ++i)
            solution[i] = i+1;
        solution[numCities] = 1;
    }

    float minDist, tempDist;
    int cityNum;
    for(int i=0; i<numCities-1; ++i)
    {
        minDist = calcEdgeDist(cities[solution[i]], cities[solution[i+1]]);
        cityNum = i+1;
        for(int j=i+2; j<numCities; ++j)
        {
            tempDist = calcEdgeDist(cities[solution[i]], cities[solution[j]]);
            if(tempDist < minDist)
            {
                minDist = tempDist;
                cityNum = j;
            }
        }
        swap(solution[i+1], solution[cityNum]);
    }
    solution[0] = solution[numCities] = 1; // starting and ending point is always labbel 1
    return solution;
}

int *TabuSearch::genRandomSolution(int *solution, bool initialize)
{
    // we can shuffle values from previous iterations
    if(initialize)
    {
        for(int i=0; i<numCities; ++i)
            solution[i] = i;
        solution[numCities] = 1;
    }
    // shuffle elements in array, without first and last element
    int r;
    for(int i=1; i<numCities-1; ++i)
    {
        r = rand()%(numCities - i - 1) + i+1; // rand in range [i+1, numCities)
        swap(solution[i], solution[r]);
    }

//    for(int i=0; i<numCities+1; ++i)
//        cout << solution[i] << " ";
//    cout << endl;

    return solution;
}

// calculates distance after swaping two edge pairs
float TabuSearch::calcNewDist(int *solution, int firstPair, int secondPair, float oldDist)
{
    if(firstPair+1 == secondPair-1) // for example: 1->2->3->4 changes to 1->3->2->4
        return (oldDist - calcEdgeDist(cities[solution[firstPair]], cities[solution[firstPair+1]]) -\
                calcEdgeDist(cities[solution[secondPair]], cities[solution[secondPair+1]]) +\
                calcEdgeDist(cities[solution[firstPair]], cities[solution[secondPair]]) +\
                calcEdgeDist(cities[solution[firstPair+1]], cities[solution[secondPair+1]]));
    else // 1->2->3->4->5 changes to 1->4->3->2->5
        return(oldDist - calcEdgeDist(cities[solution[firstPair]], cities[solution[firstPair+1]]) -\
                calcEdgeDist(cities[solution[firstPair+1]], cities[solution[firstPair+2]]) -\
                calcEdgeDist(cities[solution[secondPair-1]], cities[solution[secondPair]]) -\
                calcEdgeDist(cities[solution[secondPair]], cities[solution[secondPair+1]]) +\
                calcEdgeDist(cities[solution[firstPair]], cities[solution[secondPair]]) +\
                calcEdgeDist(cities[solution[secondPair]], cities[solution[firstPair+2]]) +\
                calcEdgeDist(cities[solution[secondPair-1]], cities[solution[firstPair+1]]) +\
                calcEdgeDist(cities[solution[firstPair+1]], cities[solution[secondPair+1]]));
}

void TabuSearch::calcTabuIters(int *localSolution, int i, int &tabuIFirst, int &tabuISec)
{
    if( localSolution[i] < localSolution[i+1] )
    {
        tabuIFirst = localSolution[i] - 1;
        tabuISec = numCities - localSolution[i+1];
    }
    else
    {
        tabuIFirst = localSolution[i+1] - 1;
        tabuISec = numCities - localSolution[i];
    }
}

void TabuSearch::calcTabuIters2(int *localSolution, int i, int j, int &tabuIFirst, int &tabuISec)
{
    if( localSolution[i] < localSolution[j] )
    {
        tabuIFirst = localSolution[i] - 1;
        tabuISec = numCities - localSolution[j];
    }
    else
    {
        tabuIFirst = localSolution[j] - 1;
        tabuISec = numCities - localSolution[i];
    }
}

void TabuSearch::tabuSearch()
{
    int *localSolution = new int[numCities+1];
    int *bestLocalSolution = new int[numCities+1];
    float globalSolutionDist, localSolutionDist, bestLocalSolutionDist;
    float currentDist;
    float tempDist;

    int firstPair, secondPair;
    int whichTry;
    int tabuIFirst, tabuISec, tabuJFirst, tabuJSec;

    for(int tries=0; tries<numMaxTries; ++tries)
    {
        // initialize best global solution, and local solutions
        if(tries == 0)
        {
            localSolution = genGreedySolution(localSolution, true);
            for(int i=0; i<numCities+1; ++i)
                globalSolution[i] = bestLocalSolution[i] = localSolution[i];
        }
        else
        {
            localSolution = genRandomSolution(localSolution, false);
            for(int i=0; i<numCities+1; ++i)
                bestLocalSolution[i] = localSolution[i];
        }
        // calculate overall distance in newly created solution
        currentDist = bestLocalSolutionDist = calcAllDist(localSolution);
        //cout << "current dist: " << currentDist << endl;
        if(tries == 0)
            globalSolutionDist = bestLocalSolutionDist;

        for(int iter=0; iter<numIterations; ++iter)
        {
            localSolutionDist = numeric_limits<float>::max();
            firstPair = secondPair = 0;
            for(int i=0; i<numCities -2; ++i) // go through all possible two edge swaps, last pair is [numCities-3, numCities-2] and [numCities-1,numCities]
            {
                for(int j=i+2; j<numCities; ++j) // go through all ending pairs, last one is [numCities-1, numCities]
                {

                    calcTabuIters2(localSolution, i, j, tabuIFirst, tabuISec);

                    if(!tabuList[tabuIFirst][tabuISec])
                    {
                        tempDist = calcNewDist(localSolution, i, j, currentDist);

                        //find best possible solution in current iteration, it can be worse than our best solution from previous iteration,
                        // but global variable will remember best one so far in a try
                        if(tempDist < localSolutionDist)
                        {
                            localSolutionDist = tempDist;
                            firstPair = i; // [i, i+1]
                            secondPair = j; // [j, j+1]
                        }
                    }
                    else //tabu
                    {
                        //check for aspire, skip tabu condition if current result is the best one so far in current trie
                        tempDist = calcNewDist(localSolution, i, j, currentDist);
                        if(tempDist < bestLocalSolutionDist)
                        {
                            localSolutionDist = tempDist;
                            firstPair = i; // [i, i+1]
                            secondPair = j; // [j, j+1]
                            //cout << "tabu aspire!" << endl;
                        }

                    }
                    // decrease tabu duration
                    if( tabuList[tabuIFirst][tabuISec] ) --tabuList[tabuIFirst][tabuISec];
                }
            }
            //cout << "przed:  " << calcAllDist(localSolution) << endl;
            // update local solution, swap second value from first pair, with first value from second pair, ex. [a,b] [c,d] -> [a,c] [b,d]

            if(secondPair == 0)
            {
                cout << "popraw parametry, secondPair==0" << endl;
            }
            calcTabuIters2(localSolution, firstPair, secondPair, tabuIFirst, tabuISec);
            tabuList[tabuIFirst][tabuISec] = tabuListLength;


            swap(localSolution[firstPair+1], localSolution[secondPair]);


            //cout << "po    : " << calcAllDist(localSolution) << endl;
            //cout << "mustbe: " << localSolutionDist << endl;
            currentDist = localSolutionDist;

            // update best local solution if found one
            if(localSolutionDist < bestLocalSolutionDist)
            {
                bestLocalSolutionDist = localSolutionDist;
                for(int i=0; i<numCities+1; ++i)
                    bestLocalSolution[i] = localSolution[i];
            }
        }
        // update best global solution
        if(bestLocalSolutionDist < globalSolutionDist)
        {
            globalSolutionDist = bestLocalSolutionDist;
            for(int i=0; i<numCities+1; ++i)
                globalSolution[i] = bestLocalSolution[i];
            whichTry = tries;
        }
    }
    cout << "global solution: " << globalSolutionDist << endl;
    //cout << "assert: " << calcAllDist(globalSolution) << endl;
    cout << "whichTry: " << whichTry << endl;

//    for(int i=0; i<numCities-1; ++i)
//    {
//        for(int j=0; j<numCities-1-i; ++j)
//            cout << tabuList[i][j] << " ";
//        cout << endl;
//    }

    delete []localSolution;
    delete []bestLocalSolution;
}

float TabuSearch::checkAllPermutations()
{
    float minDist = numeric_limits<float>::max();
    float tmpDist;
    int count=0;
    std::sort(globalSolution+1, globalSolution + numCities);
    for(int i=0; i<numCities+1; ++i)
        cout << globalSolution[i] << " ";
    do {
        tmpDist = calcAllDist(globalSolution);
        if(tmpDist < minDist) minDist = tmpDist;
        if(count==100000) cout << "LOL" << minDist << endl;
        ++count;
    }while(next_permutation(globalSolution+1, globalSolution + numCities));
    return minDist;
}


int main()
{
    srand(time(NULL));
    string fileName = "dane.txt";
    TabuSearch tabu(10, 100);
    tabu.initFromFile(fileName);
    tabu.tabuSearch();
    //tabu.printBestSolution();
    //float wynik = tabu.checkAllPermutations();
    //cout << wynik << endl;
    return 0;
}

