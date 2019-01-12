// g++ --std=c++11 optimize.cpp -o optimize.exe
// optimize in.txt out.txt

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

struct Point
{
  float x;
  float y;
};

struct Solution
{
  Solution() {}

  Solution(const std::vector<Point*>& initialData)
  {
    arr = initialData;
    calcCost();
  }

  float calcCost() {
    cost = 0.0f;
    for (int i = 0; i < arr.size() - 1; i++) {
      cost += sqrtf(
        (arr[i]->x - arr[i + 1]->x) * (arr[i]->x - arr[i + 1]->x) +
        (arr[i]->y - arr[i + 1]->y) * (arr[i]->y - arr[i + 1]->y)
        );
    }
    return cost;
  }

  Point * operator [](int i) const 
  {
    return arr[i];
  }

  Point * & operator [](int i)
  {
    return arr[i];
  }

  std::vector<Point*> arr;
  float cost;
};

class Optimizer
{
public: 

  Optimizer(): 
    m_initialized(false),
    m_loops(600),
    m_populationSize(10000),
    m_percToKeep(0.3),
    m_totalTime(0),
    m_calcCostTime(0),
    m_sortTime(0)

  {
    m_rng.seed(std::random_device()());
  }

  ~Optimizer()
  {
    for (Solution *s: m_population)
    {
      delete s;
    }
    for (Point *p: m_inData.arr)
    {
      delete p;
    }
  }

  bool init(const std::string& inPath, const std::string& outPath)
  {
    if (m_initialized) {
      return 1;
    }

    m_inPath = inPath;
    m_outPath = outPath;

    std::string line;
    std::ifstream inFile(inPath.c_str());
    if (inFile.is_open())
    {
      while (getline(inFile, line))
      {
        std::stringstream inStream(line);
        Point *point = new Point();
        inStream >> point->x >> point->y;
        m_inData.arr.push_back(point);
      }
      m_inData.calcCost();
      inFile.close();
    }
    else
    {
      return 2;
    }
    m_numberOfPoints = m_inData.arr.size();
    m_rndPoint = std::uniform_int_distribution<std::mt19937::result_type>(0, m_numberOfPoints - 1);

    for (int i = 0; i < m_populationSize; i++)
    {
      m_population.push_back(new Solution());
    }
    
    preaprePopulation();

    m_initialized = true;
    return 0;
  }

  void optimize()
  {
    clock_t globalStart = clock();
    clock_t start;
    std::vector<float> costs(m_populationSize);
    for (int loop = 0; loop < m_loops; loop++)
    {
      // calc all coasts
      start = clock();
      int minCostId = -1;
      for (int i = 0; i < m_populationSize; i++)
      {
        costs[i] = m_population[i]->calcCost();
        if (minCostId < 0 || costs[i] < costs[minCostId])
        {
          minCostId = i;
        }
      }

      // check for end (store new best solution)
      if (m_bestSolutionCost < 0.0f || costs[minCostId] < m_bestSolutionCost)
      {
        // shallow copy
        m_bestSolution = *m_population[minCostId];
        m_bestSolutionCost = costs[minCostId];
      }
      m_calcCostTime += clock() - start;

      // slecetion
      start = clock();
      std::sort(m_population.begin(), m_population.end(), 
        [](const Solution *s1, const Solution *s2) -> bool
        {
          return s1->cost > s2->cost;
        }
      );
      m_sortTime = clock() - start;

      // cross-over

      // mutation

    }
    m_totalTime = clock() - globalStart;
  }

  float getInputCost()
  {
    if (m_initialized) 
    {
      return m_inData.cost;
    }
    return -1.0f;
  }

  long getGlobalTime() { return m_totalTime; }

private:
  
  std::string m_inPath;
  std::string m_outPath;
  Solution m_inData;
  std::vector<Solution*> m_population;

  Solution m_bestSolution;
  float m_bestSolutionCost = -1.0f;

  std::mt19937 m_rng;
  std::uniform_int_distribution<std::mt19937::result_type> m_rndPoint;

  long m_totalTime;
  long m_calcCostTime;
  long m_sortTime;

  bool m_initialized;

  int m_loops;
  int m_populationSize;
  float m_percToKeep;
  int m_numberOfPoints;

  void preaprePopulation()
  {
    for (int i = 0; i < m_populationSize; i++)
    {
      for (int j = 0; j < m_numberOfPoints; j++) {
        m_population[i]->arr.push_back(m_inData[j]);
      }
      m_population[i]->calcCost();
      for (int j = 0; j < m_numberOfPoints; j++) 
      {
        std::swap((*m_population[i])[j], (*m_population[i])[m_rndPoint(m_rng)]);
      }
    }
  }

};

bool parseArguments(int argc, char **argv, std::string& inPath, std::string& outPath)
{
  if (argc != 3) {
    return 1;
  }
  inPath = argv[1];
  outPath = argv[2];
  return 0;
}

int main(int argc, char *argv[]) {

  std::string inPath;
  std::string outPath;
  
  if (parseArguments(argc, argv, inPath, outPath)) {
    std::cout << "Run: optimize [in] [out]" << std::endl;
    return 1;
  }
  
  Optimizer *optimizer = new Optimizer();
  if (int err = optimizer->init(inPath, outPath))
  {
    switch(err)
    {
    case 1: 
      std::cout << "Mltiple initialization" << std::endl;
      break;
    case 2: 
      std::cout << "Bad input file" << std::endl;
      break;
    }
    
    return 1;
  }

  std::cout << "Initial cost: " << optimizer->getInputCost() << std::endl;
  
  optimizer->optimize();

  delete optimizer;
  std::cout << "End" << std::endl;

  return 0;
}