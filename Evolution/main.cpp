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
#include <thread>


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
    friend class EvolutionTSP;
    int cityNum;
    double x;
    double y;
public:
    City(int cityNum, double x, double y) : cityNum(cityNum), x(x), y(y) { }
    City() { }
};

class EvolutionTSP
{
    int numCities;
    vector<City*> globalSolution;
    double globalSolutionDist;
    vector<City*> cities;
    //int **tabuList;
    int maxExecutionTime;

    // EVOLUTION INVER-OVER PARAMETERS
    int numMaxTries;
    int populationSize;
    int selfInversionProbability;

    vector< vector<City*> > population;
    vector< double > populationScores;

    double greedySolutionDist;
    bool finishedWaiting = false;
public:
    EvolutionTSP(); // initialization
    ~EvolutionTSP();
    void furtherInit(istream &in);
    void init(string &file);
    //City** genGreedySolution(City **solution, bool initialize); // generates greedy entry solution
    //City** genRandomSolution(City **solution, bool initialize);
    double calcAllDist(vector<City *> cityLabel) const; // returns distance of the whole route
    double calcEdgeDist(City *a, City *b) const; // returns distance between two cities
    float calcNewDist(vector<City *> solution, int firstPair, int secondPair, double oldDist) const;
    //float checkAllPermutations();
    void printBestSolution() const;
    void startEvolutionTSP();
    //void endComputations(City **solution);
    void initPopulation();
    int getMaxExecutionTime();
    void setFinishedWaiting(bool);
    void endComputations();
};

inline void EvolutionTSP::setFinishedWaiting(bool finished)
{
    finishedWaiting = finished;
}

inline int EvolutionTSP::getMaxExecutionTime()
{
    return maxExecutionTime;
}

void EvolutionTSP::initPopulation()
{
    population.reserve(populationSize);
    populationScores.reserve(populationSize);

    for (int i = 0; i < populationSize; ++i)
    {
        vector<City *> vec;
        vec.reserve(numCities);
        std::copy(cities.begin(), cities.end() - 1, std::back_inserter(vec));
        std::random_shuffle(vec.begin(), vec.end());
        populationScores.push_back(calcAllDist(vec));
        population.push_back(vec);
    }
    cout << "size" << populationScores.size() << endl;



    /*cout << "PRZED" << endl;
    for (auto iter = populationScores.begin(); iter != populationScores.end(); ++iter)
    {
        cout << "score " << *iter << endl;
    }*/

    double bestDist = populationScores[0];
    for (auto iter = populationScores.begin() + 1; iter != populationScores.end(); ++iter)
    {
        //cout << "wow " << bestDist << endl;
        bestDist = std::min(bestDist, *iter);
    }
    cout << "best score before computations: " << bestDist << endl;

    //cout << "best dist found: " << bestDist << endl;
}


/*
inline void EvolutionTSP::endComputations(City **solution)
{
    double score = calcAllDist(solution);
    if (globalSolutionDist < score) // globalSolutionDist is greedyDistance
    {
        cout << globalSolutionDist;
        //  for(int i=0; i<numCities+1; ++i)
        //      clog << globalSolution[i]->cityNum << " ";
    }
    else
    {
        cout << score;
        //  for(int i=0; i<numCities+1; ++i)
        //     clog << solution[i]->cityNum << " ";
    }
    cout << "END COMPUTATIONS" << endl;

    exit(0);
}
*/

inline void EvolutionTSP::printBestSolution() const
{
    cout << "PRINT BEST SOL" << endl;


    double bestDist = populationScores[0];
    cout << "pop: " << populationScores[0] << endl;
    int i = 1;
    int bestIter = 0;
    for (auto iter = populationScores.begin() + 1; iter != populationScores.end(); ++iter)
    {
        //cout << "wow " << bestDist << endl;
        //bestDist = std::min(bestDist, *iter);
        cout << "pop " << i <<", dist: " << *iter << endl;
        if (bestDist >= *iter)
        {
            cout << "best iter: " << i << endl;
            bestDist = *iter;
            bestIter = i;
        }
        ++i;
    }
    cout << "best iter : " << bestIter << endl;
    cout << "best dist found : " << bestDist << endl;
    cout << "assert best dist: " << calcAllDist(population[bestIter]) << endl;

    cout << "best route found: " << endl;
    for (auto && i : population[bestIter])
        cout << i->cityNum << " ";
    cout << (*population[bestIter].begin())->cityNum << endl;
}

inline void EvolutionTSP::endComputations()
{
     double bestDist = populationScores[0];
     int i = 1;
      int bestIter = 0;
     for (auto iter = populationScores.begin() + 1; iter != populationScores.end(); ++iter)
     {
         if (bestDist >= *iter)
         {
             bestDist = *iter;
             bestIter = i;
         }
         ++i;
     }
     cout << "best dist found : " << bestDist << endl;
     cout << "assert best dist: " << calcAllDist(population[bestIter]) << endl;
     exit(0);
     // TODO sprawdzanie z greedy + wypisywanie n+1 elementowego ciągu od miasta nr 1.


    /*double score = calcAllDist(solution);
    if(globalSolutionDist < score) // globalSolutionDist is greedyDistance
    {
        cout << globalSolutionDist;
       //  for(int i=0; i<numCities+1; ++i)
       //      clog << globalSolution[i]->cityNum << " ";
    }
    else
    {
        cout << score;
       //  for(int i=0; i<numCities+1; ++i)
        //     clog << solution[i]->cityNum << " ";
    }
    cout << "END COMPUTATIONS" << endl;

    exit(0); */
}

inline EvolutionTSP::EvolutionTSP()
{
    //this->numMaxTries = numMaxTries;
}

EvolutionTSP::~EvolutionTSP()
{
    // or just let it go out of scope
    vector<City*>().swap(cities); //That will create an empty vector with no memory allocated and swap it with cities, effectively deallocating the memory

}

void EvolutionTSP::init(string &fileName)
{
    if (!fileName.empty())
    {
        ifstream file(fileName.c_str());
        if (file.good())
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

void EvolutionTSP::furtherInit(istream &file)
{
    file >> numCities;

    // INVER-OVER PARAMETERS
    selfInversionProbability = 2; // 2 percent
    populationSize = 20;
    numMaxTries = 200000;


    //globalSolution = new City*[numCities + 1];
    //cities = new City*[numCities + 1];

    globalSolution.reserve(numCities + 1);
    cities.reserve(numCities + 1);

    int cityNum;
    double x, y;
    for (int i = 0; i<numCities; ++i)
    {
        file >> cityNum >> x >> y;
        cities.push_back( new City(cityNum, x, y) );
    }
    cities.push_back( cities[0] );
    file >> maxExecutionTime;
    maxExecutionTime = float(maxExecutionTime) * 0.9;
}

inline double EvolutionTSP::calcEdgeDist(City *a, City *b) const
{
    return sqrt(pow((a->x - b->x), 2) + pow((a->y - b->y), 2));
}

inline double EvolutionTSP::calcAllDist(vector<City*> cityLabel) const
{
    double dist = 0.0;
    //for (int i = 0; i<numCities; ++i) // go through all edges, city's length is numCities+1
    //{
    //	dist += calcEdgeDist(cityLabel[i], cityLabel[i + 1]);
    //}
    for(auto it = cityLabel.begin(); it != cityLabel.end() - 2; ++it) // TODO poprawic to, bo retarded sie wydaje
        dist += calcEdgeDist( *it, *(it + 1) );
    dist += calcEdgeDist(*(cityLabel.end() - 1), *(cityLabel.begin()) ); // connection from last city to first one
    return dist;
}

/*
// generate greedy solution based on Cities **cities
City** EvolutionTSP::genGreedySolution(City **solution, bool initialize)
{
    if (initialize)
    {
        for (int i = 0; i<numCities + 1; ++i)
            solution[i] = cities[i];
    }

    float minDist, tempDist;
    int cityNum;
    for (int i = 0; i<numCities - 1; ++i)
    {
        minDist = calcEdgeDist(solution[i], solution[i + 1]);
        cityNum = i + 1;
        for (int j = i + 2; j<numCities; ++j)
        {
            tempDist = calcEdgeDist(solution[i], solution[j]);
            if (tempDist < minDist)
            {
                minDist = tempDist;
                cityNum = j;
            }
        }
        swap(solution[i + 1], solution[cityNum]);
    }
    solution[0] = solution[numCities] = cities[0]; // starting and ending point is always labbel 1
    return solution;
}

City **EvolutionTSP::genRandomSolution(City **solution, bool initialize)
{
    // we can shuffle values from previous iterations
    if (initialize)
    {
        for (int i = 0; i<numCities + 1; ++i)
            solution[i] = cities[i];
    }
    // shuffle elements in array, without first and last element
    int r;
    for (int i = 1; i<numCities - 1; ++i)
    {
        r = rand() % (numCities - i - 1) + i + 1; // rand in range [i+1, numCities)
        swap(solution[i], solution[r]);
    }

    //    for(int i=0; i<numCities+1; ++i)
    //        cout << solution[i] << " ";
    //    cout << endl;
    return solution;
}
*/

// calculates distance after swaping two edge pairs
inline float EvolutionTSP::calcNewDist(vector<City*> solution, int firstPair, int secondPair, double oldDist) const
{
    --firstPair; // changed from previous solution when I passed [i, i+1] and firstPair was equal i, now I pass i+1, kek, i dont want to change all other stuff for now
    if (firstPair + 1 == secondPair - 1) // for example: 1->2->3->4 changes to 1->3->2->4
        return (oldDist - calcEdgeDist(solution[firstPair], solution[firstPair + 1]) - \
            calcEdgeDist(solution[secondPair], solution[secondPair + 1]) + \
            calcEdgeDist(solution[firstPair], solution[secondPair]) + \
            calcEdgeDist(solution[firstPair + 1], solution[secondPair + 1]));
    else // 1->2->3->4->5 changes to 1->4->3->2->5
        return(oldDist - calcEdgeDist(solution[firstPair], solution[firstPair + 1]) - \
            calcEdgeDist(solution[firstPair + 1], solution[firstPair + 2]) - \
            calcEdgeDist(solution[secondPair - 1], solution[secondPair]) - \
            calcEdgeDist(solution[secondPair], solution[secondPair + 1]) + \
            calcEdgeDist(solution[firstPair], solution[secondPair]) + \
            calcEdgeDist(solution[secondPair], solution[firstPair + 2]) + \
            calcEdgeDist(solution[secondPair - 1], solution[firstPair + 1]) + \
            calcEdgeDist(solution[firstPair + 1], solution[secondPair + 1]));
}

void threadTimeChecker(EvolutionTSP *evol)
{
    int waitTime = evol->getMaxExecutionTime() * 0.92;
    std::chrono::seconds duration(waitTime);
    std::this_thread::sleep_for(duration);
    evol->setFinishedWaiting(true);
}

inline bool isCitiesNextToEachOther(vector<City*>::iterator firstCityIter, vector<City*>::iterator secondCityIter, vector<City*> &firstMember)
{
    if (firstCityIter == firstMember.begin()) // at beginning of cycle
        return (firstCityIter + 1 == secondCityIter || firstMember.end() - 1 == secondCityIter);
    else if (firstCityIter == firstMember.end() - 1) // at end of cycle
        return (firstMember.begin() == secondCityIter || firstCityIter - 1 == secondCityIter);
    else // in the middle of cycle
        return (firstCityIter + 1 == secondCityIter || firstCityIter - 1 == secondCityIter);
}


/* function reverses vector between two iterators in a cycle manner, in [upper; lower],
    example 1 - 2 - 3 - 4 - 5 - 6, reverse from 4 to 1 yields 4 - 2 - 3 - 1 - 6 - 5
*/
void inverseCities(vector<City*>::iterator lower, vector<City*>::iterator upper, vector<City*> &vec)
{
    int lowerCountElements = lower - vec.begin() + 1;
    int upperCountElements = vec.end() - upper;
    bool lowerBound = false, upperBound = false;

    if (lowerCountElements < upperCountElements)
    {
        lowerBound = true;
    }
    else if (lowerCountElements > upperCountElements)
    {
        upperBound = true;
    }
    // else we're done after first loop
    int firstReverseCount = std::min(lowerCountElements, upperCountElements);

    while (firstReverseCount > 1)
    {
        //if (upper == vec.end()) cout << "------------" << firstReverseCount << endl;
        //cout << "++++++++" << firstReverseCount << endl;

        //cout << "LOWER COUNT: " << lowerCountElements << endl;
        //cout << "UPPER COUNT: " << upperCountElements << endl;
        //if (lower == vec.end()) cout << "LOWER AT END" << endl;

        swap(*lower, *upper);
        //cout << "after swap" << endl;
        if(lower != vec.begin())
            --lower;
        //else cout << "LOWER SPADA" << endl;
        //if (upper == vec.end() - 1 && firstReverseCount >= 2) cout << "UPPER: revCount: " << firstReverseCount << endl;
        ++upper;

        --firstReverseCount;
    }
    //swap(*lower, *upper);


    if (lowerBound) // lower iterator moves to the end of vector
    {
        //cout << "NOT1" << endl;
        lower = vec.end();
        std::reverse(upper, lower);
    }
    else if(upperBound) // upper iterator moves to beginning of vector
    {
        //cout << "NOT2" << endl;
        upper = vec.begin();
        ++lower; // because reverse operates on [upper, lower) interval
        std::reverse(upper, lower);
    }
}


void EvolutionTSP::startEvolutionTSP()
{
    initPopulation();
    vector<City*> firstMember;
    firstMember.resize(numCities);
    vector<City*> *secondMember;
    vector<City*>::iterator secondCityIter, firstCityIter;
    int firstRandCity, secondRandCity;
    int repeatTimes;
    int populationMemberID;

    for (int iter = 0; iter < numMaxTries; ++iter)
    {
        //cout << "Try: " << iter << endl;
        populationMemberID = 0;
        for (auto && actualPopulationMember : population) // for each population member
        {
            std::copy(actualPopulationMember.begin(), actualPopulationMember.end(), firstMember.begin() );
            firstRandCity = rand() % numCities; // [0, numCities)
            if(finishedWaiting) endComputations();
            repeatTimes = 0;
            do
            {
                //cout << "repeatTimes: " << repeatTimes << endl;
                ++repeatTimes;
                if (rand() % 101 <= selfInversionProbability) // do self inverse
                {
                    // get different second city
                    do
                    {
                        secondRandCity = rand() % numCities;
                    } while (firstRandCity == secondRandCity);
                    firstCityIter = firstMember.begin() + firstRandCity;
                    secondCityIter = firstMember.begin() + secondRandCity;

                    if(isCitiesNextToEachOther(firstCityIter, secondCityIter, firstMember))
                        break;

                        // reverse cities between [firstIterator+1 and secondIterator]
                    if (firstCityIter < secondCityIter)
                    {
                        secondRandCity = firstCityIter + 1 - firstMember.begin();
                        std::reverse(firstCityIter + 1, secondCityIter + 1);
                    }
                    else
                    {
                        if (firstCityIter == firstMember.end() - 1)
                            secondRandCity = 0; // beginning of cycle
                        else // it's after our iterator
                            secondRandCity = firstCityIter + 1 - firstMember.begin();

                        inverseCities(secondCityIter + 1, firstCityIter, firstMember);
                    }
                }
                else // inverse range determined from another population member
                {
                    // get different population member
                    do
                    {
                        secondMember = & population[rand() % populationSize];
                    } while ( (*secondMember) == actualPopulationMember);

                    // secondCity now becomes the city "after" firstCity in secondTempMember
                    if ((secondCityIter = std::find( (*secondMember).begin(), (*secondMember).end(), *(firstMember.begin() + firstRandCity) )) != (*secondMember).end())
                    {
                        //city found
                        //cout << "CURRENT CITY SHOULD BE: " << actualPopulationMember[firstRandCity]->cityNum << endl;
                        //cout << "CURRENT CITY FOUND    : " << (*secondCityIter)->cityNum << endl;
                        if (secondCityIter == (*secondMember).end() - 1) // following city is at cycle start
                        {
                            secondCityIter = (*secondMember).begin();
                        }
                        else
                        {
                            ++secondCityIter;
                        }
                        //cout << "NEXT CITY IS      : " << (*secondCityIter)->cityNum << endl;

                        // find secondCity in our first member
                        if ((secondCityIter = std::find(firstMember.begin(), firstMember.end(), *secondCityIter)) != firstMember.end())
                        {
                            //	cout << "ASSERT NEXT CITY : " << (*secondCityIter)->cityNum << endl;
                            // city found
                            firstCityIter = firstMember.begin() + firstRandCity;
                            //	cout << "FIRST CITY ITER IS: " << (*firstCityIter)->cityNum << endl;
                            //cout << "FIRST CITY ITER SHOULD BE: " << actualPopulationMember[firstRandCity]->cityNum << endl;
                            // if secondCity is previous or next from firstCity, leave loop
                            if (isCitiesNextToEachOther(firstCityIter, secondCityIter, firstMember))
                                break;

                            // reverse cities between [firstIterator+1 and secondIterator]
                            if (firstCityIter < secondCityIter)
                            {
                                secondRandCity = firstCityIter + 1 - firstMember.begin();
                                std::reverse(firstCityIter + 1, secondCityIter + 1);
                            }
                            else
                            {
                                if (firstCityIter == firstMember.end() - 1)
                                    secondRandCity = 0; // beginning of cycle
                                else // it's after our iterator
                                    secondRandCity = firstCityIter + 1 - firstMember.begin();

                                inverseCities(secondCityIter + 1, firstCityIter, firstMember);
                            }
                        }
                        else
                        {
                            cout << "nie znaleziono miasta 2" << endl;
                        }
                    }
                    else
                    {
                        // city not found
                        cout << "nie znaleziono miasta 1" << endl;
                    }
                }
                firstRandCity = secondRandCity;

            } while(repeatTimes <= 200);

            //cout << "all dist:" << calcAllDist(firstMember) << endl;
            double dist = calcAllDist(firstMember);
            if (dist <= populationScores[populationMemberID])
            {
                std::copy(firstMember.begin(), firstMember.end(), actualPopulationMember.begin());
                populationScores[populationMemberID] = dist;
            }

            ++populationMemberID;
            //firstMember.clear();
        }
    }

    //for (auto iter = populationScores.begin(); iter != populationScores.end(); ++iter)
    //{
    //	cout << "score " << *iter << endl;
    //}


    cout << "PO " << endl;
    double bestDist = populationScores[0];
    for (auto iter = populationScores.begin() + 1; iter != populationScores.end(); ++iter)
    {
        //cout << "wow " << bestDist << endl;
        bestDist = std::min(bestDist, *iter);
    }
    cout << "best dist found: " << bestDist << endl;

}


int main()
{
    cout << "start" << endl;
    startTime = clock();
    srand(time(NULL));

    string fileName = "mar_506.txt"; //TODO WYWALIC STRINGA NA NULLA XD
    EvolutionTSP *sim = new EvolutionTSP();
    sim->init(fileName);
    std::thread checker(threadTimeChecker, sim);
    //cout << "wczytalem" << endl;
    sim->startEvolutionTSP();
    //sim.printBestSolution();
    //cout << "\n\n\n\n\n";

    cout << "end" << endl;
  //  checker.join();
    exit(0);
   // return 0;
}
