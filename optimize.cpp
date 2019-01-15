// g++ -O2 --std=c++11 optimize.cpp -o optimize.exe
// gdb32 --args optimize.exe in.txt out.txt
// optimize in.txt out.txt

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

struct Point
{
  float x;
  float y;
  int id;
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
    if (arr.size() == 0) 
    {
      return 0.0f;
    }
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

  Optimizer(bool debug): 
    m_initialized(false),
    m_loops(600),
    m_populationSize(6000),
    m_debug(debug),
    m_percToKeep(0.3),
    m_percMutation(0.06),
    m_numberOfMutations(11),
    m_totalTime(0),
    m_calcCostTime(0),
    m_sortTime(0),
    m_crossOverTime(0),
    m_populationOld(new std::vector<Solution*>()),
    m_populationNew(new std::vector<Solution*>())
  {
    m_rng.seed(time(nullptr));
    m_rndSolution = std::uniform_int_distribution<std::mt19937::result_type>(0, m_populationSize * m_percToKeep - 1);
    m_rndMutation = std::uniform_real_distribution<float>(0.0f, 1.0f);
  }

  ~Optimizer()
  {
    for (Solution *s: *m_populationOld)
    {
      delete s;
    }
    for (Solution *s: *m_populationNew)
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
      int ind = 0;
      while (getline(inFile, line))
      {
        std::stringstream inStream(line);
        Point *point = new Point();
        point->id = ind++;
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
    m_bestSolution = m_inData;
    m_rndPoint = std::uniform_int_distribution<std::mt19937::result_type>(0, m_numberOfPoints - 1);

    for (int i = 0; i < m_populationSize; i++)
    {
      m_populationOld->push_back(new Solution());
      m_populationNew->push_back(new Solution());
    }

    for (int i = 0; i < m_numberOfPoints; i++)
    {
      m_visited.push_back(false);
    }
    
    preaprePopulation();

    m_initialized = true;
    return 0;
  }

  void optimize()
  {
    globalStart = clock();
    clock_t start;
    std::vector<float> costs(m_populationSize);
    for (int loop = 0; loop < m_loops; loop++)
    {
      if (loop % (m_loops * 5 / 100) == 0)
      {
        std::cout << loop * 100 / m_loops << "% \n";
        if (m_debug)
        {
          for (int i = 0; i < 20; i++)
          {
            std::cout << (*m_populationOld)[i]->cost << " ";
          }
          for (int i = 20; i > 0; i--)
          {
            std::cout << (*m_populationOld)[m_percToKeep * m_populationSize - i]->cost << " ";
          }
        }
      }
      // calc all coasts
      start = clock();
      int minCostId = -1;
      for (int i = 0; i < m_populationSize; i++)
      {
        costs[i] = (*m_populationOld)[i]->calcCost();
        if (minCostId < 0 || costs[i] < costs[minCostId])
        {
          minCostId = i;
        }
      }

      // check for end (store new best solution)
      if (m_bestSolution.cost < 0.0f || costs[minCostId] < m_bestSolution.cost)
      {
        // shallow copy
        m_bestSolution = *(*m_populationOld)[minCostId];
        dump();
        std::cout << costs[minCostId] << std::endl;
      }
      m_calcCostTime += clock() - start;

      // slecetion
      start = clock();
      std::sort(m_populationOld->begin(), m_populationOld->end(), 
        [](const Solution *s1, const Solution *s2) -> bool
        {
          return s1->cost < s2->cost;
        }
      );
      m_sortTime += clock() - start;

      // cross-over
      start = clock();
      for (int i = 0; i < m_populationSize; i++)
      {
        Solution *solution = (*m_populationNew)[i];

        Solution *s[] = { 
          (*m_populationOld)[m_rndSolution(m_rng)], 
          (*m_populationOld)[m_rndSolution(m_rng)]
        };
        int sInd[] = { 0, 0 };

        int ind = 0;
        int sel = 0;
        //clearVisited();
        bool *visited = new bool[m_numberOfPoints]{ false };
        int streak = 0;
        float crossOvewScale = 0.1f + m_rndMutation(m_rng) * 0.7f;
        int streakMin = 1 + m_rndMutation(m_rng) * 2;
        while (ind < m_numberOfPoints)
        {
          if (streak > streakMin) 
          {
            if (m_rndMutation(m_rng) < crossOvewScale) { sel = (sel + 1) % 2; streak = 0; }
            if (sInd[sel] == m_numberOfPoints) { sel = (sel + 1) % 2; streak = 0; }
          }
          streak++;
          if (!visited[s[sel]->arr[sInd[sel]]->id])
          {
            solution->arr[ind++] = s[sel]->arr[sInd[sel]];
            visited[s[sel]->arr[sInd[sel]]->id] = true;
          }
          sInd[sel]++;
        }
        delete [] visited;
        if (int id = checkValid(solution))
        {
          std::cout << "cross-over error " << i << std::endl;
          for (int j = 0; j < m_numberOfPoints; j++)
          {
            if (j == id)
            {
              std::cout << "error ";
            }
            std::cout << (*solution)[j]->x << " " << (*solution)[j]->y << std::endl;
          }
          exit(1);
        }
      }
      m_crossOverTime += clock() - start;

      // mutation
      start = clock();
      for (int i = 0; i < m_populationSize; i++)
      {
        Solution *s = (*m_populationNew)[i];
        int mutations = m_rndMutation(m_rng) * m_numberOfMutations;
        for (int i = 0; i < mutations; i++)
        {
          int id1 = m_rndPoint(m_rng);
          int id2 = (id1 + 1) % m_numberOfPoints;
          std::swap(s->arr[id1], s->arr[id2]);
        }
        if (int id = checkValid(s))
        {
          std::cout << "mutation error " << std::endl;
          for (int j = 0; j < m_numberOfPoints; j++)
          {
            if (j == id)
            {
              std::cout << "error ";
            }
            std::cout << (*s)[j]->x << " " << (*s)[j]->y << std::endl;
          }
          exit(1);
        }
      }
      m_mutationTime += clock() - start;

      *(*m_populationNew)[0] = *(*m_populationOld)[0];

      std::swap(m_populationOld, m_populationNew);
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
  long getCalcCostTime() { return m_calcCostTime; }
  long getSortTime() { return m_sortTime; }
  long getCrossOverTime() { return m_crossOverTime; }
  long getMutationTime() { return m_mutationTime; }

  Solution &getBestSolution() { return m_bestSolution; }

private:
  
  std::string m_inPath;
  std::string m_outPath;
  Solution m_inData;
  std::vector<Solution*>* m_populationOld;
  std::vector<Solution*>* m_populationNew;
  std::vector<bool> m_visited;

  Solution m_bestSolution;

  std::mt19937 m_rng;
  std::uniform_int_distribution<std::mt19937::result_type> m_rndPoint;
  std::uniform_int_distribution<std::mt19937::result_type> m_rndSolution;
  std::uniform_real_distribution<float> m_rndMutation;

  clock_t globalStart;
  long m_totalTime;
  long m_calcCostTime;
  long m_sortTime;
  long m_crossOverTime;
  long m_mutationTime;

  bool m_initialized;

  int m_loops;
  int m_populationSize;
  float m_percToKeep;
  float m_percMutation;
  int m_numberOfMutations;
  bool m_debug;
  int m_numberOfPoints;

  void preaprePopulation()
  {
    for (int i = 0; i < m_populationSize; i++)
    {
      for (int j = 0; j < m_numberOfPoints; j++) {
        (*m_populationOld)[i]->arr.push_back(m_inData[j]);
        (*m_populationNew)[i]->arr.push_back(m_inData[j]);
      }
      // if (i < 1)
      // {
      //   continue;
      // }
      for (int j = 0; j < m_numberOfPoints; j++) 
      {
        std::swap((*(*m_populationOld)[i])[j], (*(*m_populationOld)[i])[m_rndPoint(m_rng)]);
      }
      if (int id = checkValid((*m_populationOld)[i]))
      {
        std::cout << "preparation error " << i << std::endl;
        for (int j = 0; j < m_numberOfPoints; j++)
        {
          if (j == id)
          {
            std::cout << "error ";
          }
          std::cout << (*(*m_populationOld)[i])[j]->x << " " << (*(*m_populationOld)[i])[j]->y << std::endl;
        }
        exit(1);
      }
    }
  }

  int checkValid(Solution* s)
  {
    if (!m_debug) 
    {
      return 0;
    }
    std::unordered_set<Point*> used;
    for (int i = 0; i < m_numberOfPoints; i++) 
    {
      if (used.find(s->arr[i]) != used.end())
      {
        return i;
      }
      used.insert(s->arr[i]);
    }
    return 0;
  }

  bool dump()
  {
    std::ofstream outFile(m_outPath.c_str());
    if (!outFile.is_open())
    {
      return false;
    }
    outFile << m_bestSolution.cost << std::endl;
    outFile << 1.0f * (clock() - globalStart) / CLOCKS_PER_SEC << std::endl;
    for (Point *p: m_bestSolution.arr)
    {
      outFile << p->x << " " << p->y << std::endl;
    }
    outFile.close();
    return true;
  }

  void clearVisited() {
    for (int i = 0; i < m_numberOfPoints; i++)
    {
      m_visited[i] = false;
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
  
  Optimizer *optimizer = new Optimizer(false);
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

  std::cout << "Accumulated time: " << optimizer->getGlobalTime() / 1000.0f << std::endl;
  std::cout << "Calc cost time: " << optimizer->getCalcCostTime() / 1000.0f << std::endl;
  std::cout << "Sorting time: " << optimizer->getSortTime() / 1000.0f << std::endl;
  std::cout << "Cross-over time: " << optimizer->getCrossOverTime() / 1000.0f << std::endl;
  std::cout << "Mutation time: " << optimizer->getMutationTime() / 1000.0f << std::endl;

  Solution &bestSolution = optimizer->getBestSolution();
  std::cout << "Best solution cost: " << bestSolution.cost << std::endl;

  delete optimizer;
  std::cout << "End" << std::endl;

  return 0;
}