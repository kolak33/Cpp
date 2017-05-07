#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>   // rand
#include <math.h>    // sqrt, pow, constant e
#include <limits>    // numeric_limits<float>::max()
#include <ctime>
#include <algorithm>
#include <vector>
#include <utility>
using std::string; using std::cin;
using std::cout; using std::endl;
using std::ifstream; using std::istream;
using std::clog; using std::swap;
using std::numeric_limits; using std::next_permutation;
using std::vector; using std::pair;

clock_t startTime;
clock_t actualTime;

class City
{
private:
    friend class SimulatedAnnealing;
    int cityNum;
    double x;
    double y;
public:
    City(int cityNum, double x, double y) : cityNum(cityNum), x(x), y(y) { }
    City() { }
};

class SimulatedAnnealing
{
private:
    int simListLength;
    int numCities;
    int numMaxTries;
    City **globalSolution;
    double globalSolutionDist;
    City **cities;
    int **tabuList;
    bool initialized;

    double endTemperatureBound;
    double temperature;
    double annealingFactor;
    int maxExecutionTime;
    double greedySolutionDist;
public:
    SimulatedAnnealing(int numMaxTries); // initialization
    ~SimulatedAnnealing();
    void furtherInit(istream &in);
    void init(string &file);
    City** genGreedySolution(City **solution, bool initialize); // generates greedy entry solution
    City** genRandomSolution(City **solution, bool initialize);
    double calcAllDist(City **cityLabel) const; // returns distance of the whole route
    double calcEdgeDist(City *a, City *b) const; // returns distance between two cities
    float calcNewDist( City **solution, int firstPair, int secondPair, double oldDist) const;
    void calcTabuIters(City **localSolution, int i, int &tabuIFirst, int &tabuISec);
    void calcTabuIters2(City **localSolution, int i, int j, int &tabuIFirst, int &tabuISec);
    float checkAllPermutations();
    void printBestSolution() const;
    void simmulatedAnnealing();
    void endComputations(City **solution);
    bool boltzmanCondition(double temperature, double oldDistance, double newDistance);
};

inline bool SimulatedAnnealing::boltzmanCondition(double temperature, double oldDistance, double newDistance)
{
    return (double(rand()) / double(RAND_MAX) < exp(( -newDistance + oldDistance ) / temperature));
}

inline void SimulatedAnnealing::endComputations(City **solution)
{
    double score = calcAllDist(solution);
    if(globalSolutionDist < score) // globalSolutionDist is greedyDistance
    {
        cout << globalSolutionDist;
        // for(int i=0; i<numCities+1; ++i)
        //     clog << globalSolution[i]->cityNum << " ";
    }
    else
    {
        cout << score;
        // for(int i=0; i<numCities+1; ++i)
        //     clog << solution[i]->cityNum << " ";
    }
    cout << "END COMPUTATIONS" << endl;

    exit(0);
}

inline void SimulatedAnnealing::printBestSolution() const
{
    for(int i=0; i<numCities+1; ++i)
        clog << globalSolution[i]->cityNum << " ";
}

inline SimulatedAnnealing::SimulatedAnnealing(int numMaxTries)
{
    this->numMaxTries = numMaxTries;
    initialized = false;
}

SimulatedAnnealing::~SimulatedAnnealing()
{
    if(initialized)
    {
        for(int i=0; i<numCities-1; ++i)
        {
            delete cities[i];
            delete tabuList[i];
        }
        delete cities[numCities-1]; delete cities[numCities];
        delete []globalSolution;
        delete []tabuList;
        delete []cities;
    }

}

void SimulatedAnnealing::init(string &fileName)
{
    if(!fileName.empty())
    {
        ifstream file(fileName.c_str());
        if(file.good())
        {
            furtherInit(file);
            file.close();
        }
        else
            clog << "error opening file" << endl;
    }
    else
       furtherInit(cin);
}

void SimulatedAnnealing::furtherInit(istream &file)
{
        file >> numCities;


        if(numCities >= 6000)
        {
            temperature = 1;
            annealingFactor = 0.9999;
        }
        else
        {
            temperature = 1;
            annealingFactor = 0.99999;
        }
        endTemperatureBound = 0.00001;

        globalSolution = new City*[numCities+1];
        cities = new City*[numCities+1];
        tabuList = new int *[numCities-1];
        for(int i=0; i<numCities-1; ++i)
        {
            tabuList[i] = new int[numCities-1-i];
            for(int j=0; j<numCities-1-i; ++j)
                tabuList[i][j] = 0;
        }

        int cityNum;
        double x, y;
        for(int i=0; i<numCities; ++i)
        {
           file >> cityNum >> x >> y;
           cities[i] = new City(cityNum, x, y);
        }
        cities[numCities] = cities[0];
        file >> maxExecutionTime;
        maxExecutionTime = float(maxExecutionTime) * 0.9;
        initialized = true;
}

inline double SimulatedAnnealing::calcEdgeDist(City *a, City *b) const
{
    return sqrt(pow((a->x - b->x), 2) + pow((a->y - b->y), 2));
}

inline double SimulatedAnnealing::calcAllDist(City **cityLabel) const
{
    double dist = 0.0;
    for(int i=0; i<numCities; ++i) // go through all edges, city's length is numCities+1
    {
        dist += calcEdgeDist(cityLabel[i], cityLabel[i+1]);
    }
    return dist;
}

// generate greedy solution based on Cities **cities
City** SimulatedAnnealing::genGreedySolution(City **solution, bool initialize)
{
    if(initialize)
    {
        for(int i=0; i<numCities+1; ++i)
            solution[i] = cities[i];
    }

    float minDist, tempDist;
    int cityNum;
    for(int i=0; i<numCities-1; ++i)
    {
        minDist = calcEdgeDist(solution[i], solution[i+1]);
        cityNum = i+1;
        for(int j=i+2; j<numCities; ++j)
        {
            tempDist = calcEdgeDist(solution[i], solution[j]);
            if(tempDist < minDist)
            {
                minDist = tempDist;
                cityNum = j;
            }
        }
        swap(solution[i+1], solution[cityNum]);
    }
    solution[0] = solution[numCities] = cities[0]; // starting and ending point is always labbel 1
    return solution;
}

City **SimulatedAnnealing::genRandomSolution(City **solution, bool initialize)
{
    // we can shuffle values from previous iterations
    if(initialize)
    {
        for(int i=0; i<numCities+1; ++i)
            solution[i] = cities[i];
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
inline float SimulatedAnnealing::calcNewDist(City **solution, int firstPair, int secondPair, double oldDist) const
{
    --firstPair; // changed from previous solution when I passed [i, i+1] and firstPair was equal i, now I pass i+1, kek, i dont want to change all other stuff for now
    if(firstPair+1 == secondPair-1) // for example: 1->2->3->4 changes to 1->3->2->4
        return (oldDist - calcEdgeDist(solution[firstPair], solution[firstPair+1]) -\
                calcEdgeDist(solution[secondPair], solution[secondPair+1]) +\
                calcEdgeDist(solution[firstPair], solution[secondPair]) +\
                calcEdgeDist(solution[firstPair+1], solution[secondPair+1]));
    else // 1->2->3->4->5 changes to 1->4->3->2->5
        return(oldDist - calcEdgeDist(solution[firstPair], solution[firstPair+1]) -\
                calcEdgeDist(solution[firstPair+1], solution[firstPair+2]) -\
                calcEdgeDist(solution[secondPair-1], solution[secondPair]) -\
                calcEdgeDist(solution[secondPair], solution[secondPair+1]) +\
                calcEdgeDist(solution[firstPair], solution[secondPair]) +\
                calcEdgeDist(solution[secondPair], solution[firstPair+2]) +\
                calcEdgeDist(solution[secondPair-1], solution[firstPair+1]) +\
                calcEdgeDist(solution[firstPair+1], solution[secondPair+1]));
}

inline void SimulatedAnnealing::calcTabuIters(City **localSolution, int i, int &tabuIFirst, int &tabuISec)
{
    if( localSolution[i]->cityNum < localSolution[i+1]->cityNum )
    {
        tabuIFirst = localSolution[i]->cityNum - 1;
        tabuISec = numCities - localSolution[i+1]->cityNum;
    }
    else
    {
        tabuIFirst = localSolution[i+1]->cityNum - 1;
        tabuISec = numCities - localSolution[i]->cityNum;
    }
}

inline void SimulatedAnnealing::calcTabuIters2(City **localSolution, int i, int j, int &tabuIFirst, int &tabuISec)
{
    if( localSolution[i]->cityNum < localSolution[j]->cityNum )
    {
        tabuIFirst = localSolution[i]->cityNum - 1;
        tabuISec = numCities - localSolution[j]->cityNum;
    }
    else
    {
        tabuIFirst = localSolution[j]->cityNum - 1;
        tabuISec = numCities - localSolution[i]->cityNum;
    }
}


void SimulatedAnnealing::simmulatedAnnealing()
{
    City **localSolution = new City*[numCities+1];
    double localSolutionDist;

    double tempDist;
    double currentTemperature;

    for(int tries=0; tries<numMaxTries; ++tries)
    {
        currentTemperature = temperature;
        if(tries == 0)
        {
            localSolution = genGreedySolution(localSolution, true);
            for(int i=0; i<numCities+1; ++i)
                globalSolution[i] = localSolution[i];
            cout << "GREEDY: " << calcAllDist(globalSolution) << endl;
            localSolutionDist = globalSolutionDist = calcAllDist(localSolution);
        }
        else
        {
            localSolution = genRandomSolution(localSolution, false);
            localSolutionDist = calcAllDist(localSolution);
            cout << "STARTING SOLUTION: " << calcAllDist(localSolution) << endl;
        }


        while(currentTemperature > endTemperatureBound) // end computations condition
        { 
            for(int iterAtTemp = 0; iterAtTemp < 100; ++iterAtTemp)
            {
                int firstCity = rand()%(numCities-1) + 1; // [1; numCities-1]
                int secondCity = rand()%(numCities-1) + 1; // [1; numCities-1]
                if(firstCity > secondCity) swap(firstCity, secondCity);

                if(firstCity != secondCity)
                {
                    tempDist = calcNewDist(localSolution, firstCity, secondCity, localSolutionDist);

                    if(tempDist < localSolutionDist)
                    {
                        swap(localSolution[firstCity], localSolution[secondCity]);
                        localSolutionDist = tempDist;
                    }
                    else if( boltzmanCondition(currentTemperature, localSolutionDist, tempDist) )
                    {
                        swap(localSolution[firstCity], localSolution[secondCity]);
                        localSolutionDist = tempDist;
                    }
                }
            }

            actualTime = clock();
            if( (actualTime - startTime) / CLOCKS_PER_SEC >= maxExecutionTime) endComputations(localSolution);
            currentTemperature = currentTemperature * annealingFactor; // update temperature in each iteration
        }

        cout << "LOCAL DIST       : " << calcAllDist(localSolution) << endl;
        cout << "LOCAL DIST ASSERT: " << localSolutionDist << endl;
        cout << "\n\n";
    }
//    cout << "global solution: " << globalSolutionDist << endl;
//    cout << "assert: " << calcAllDist(globalSolution) << endl;
//    cout << "whichTry: " << whichTry << endl;

    delete []localSolution;
}

float SimulatedAnnealing::checkAllPermutations()
{
    float minDist = numeric_limits<float>::max();
    float tmpDist;
    int count=0;
    std::sort(globalSolution+1, globalSolution + numCities);

    while(next_permutation(globalSolution+1, globalSolution + numCities))
    {
        tmpDist = calcAllDist(globalSolution);
        if(tmpDist < minDist)
            minDist = tmpDist;

        if(count==20000000) {count=0; cout << "min so far: " << minDist << endl;}
        ++count;
    }

    return minDist;
}


int main()
{
    startTime = clock();
    srand(time(NULL));

    string fileName = "test-data4.txt"; //TODO WYWALIC STRINGA NA NULLA XD
    SimulatedAnnealing sim(3);
    sim.init(fileName);
    //cout << "wczytalem" << endl;
    sim.simmulatedAnnealing();
    //sim.printBestSolution();
    cout << "\n\n\n\n\n";


    return 0;
}
