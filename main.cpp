#include <iostream>
#include <iomanip> /* setprecision */
#include <chrono> /* to calculate the execution time */
#include <time.h> /* for function time(0) used for the seed */
#include <cmath> /* math functions, such as sin, pow, ceil etc*/
#include <vector> /* obviously, vector */
#include <random> /* for mt19937 */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

// Global parameters
int dimension = 30;
double minimum = NULL;
double maximum = NULL;
int precision = 10;
int requiredBits = NULL; // no. of bits required to represent a SINGLE number
vector<bool> bitstring; // will have the size of requiredBits * dimension, therefore being able to contain all required numbers for a single function evaluation
vector<double> converted; // converted numbers (from the bitstring vector candidates)
int maxIterations=10;



mt19937 randomNumber(time(0));
uniform_int_distribution<int> distribution(1, requiredBits);

enum functions
{
	DeJong = 1,
	Schwefel,
	Rastrigin,
	Michalewicz
};

long double DeJongFunction()
{
	double sum = 0;

	for (int i = 0; i < dimension; i++)
		sum += converted[i] * converted[i];

	return sum;
}

long double SchwefelFunction()
{
	double sum = 0;

	for (int i = 0; i < dimension; i++)
		sum -= converted[i] * sin(sqrt(abs(converted[i])));

	return 418.9829 * dimension + sum;
}

long double RastriginFunction()
{
	double sum = 10 * dimension;

	for (int i = 0; i < dimension; i++)
		sum += converted[i] * converted[i] - 10 * cos(2 * M_PI * converted[i]);

	return sum;
}

long double MichalewiczFunction()
{
	double sum = 0;
	int m = 10;

	for (int i = 0; i < dimension; i++)
		sum -= sin(converted[i]) * pow(sin((i + 1) * converted[i] * converted[i] / M_PI), 2 * m);

	return sum;
}

void setParameters(functions& currentFunction)
{
	// resetting the two vectors if needed
	while (!bitstring.empty())
	{
		bitstring.pop_back();
	}
	while (!converted.empty())
	{
		converted.pop_back();
	}
	switch (currentFunction)
	{
	case functions::DeJong:
		minimum = -5.12;
		maximum = 5.12;
		break;
	case functions::Schwefel:
		minimum = -600;
		maximum = 600;
		break;
	case functions::Rastrigin:
		minimum = -5.12;
		maximum = 5.12;
		break;
	case functions::Michalewicz:
		minimum = 0;
		maximum = 3.14;
		break;
	default:
		std::cout << "Error at setting the function domain";
		break;
	}
	requiredBits = std::ceil(log2((maximum - minimum) * pow(10, precision)));
}


vector<bool> generateRandBitNumber(int lungime, vector<bool> &biti) {

  for (int i = 0; i < lungime; i++) {
    biti.push_back(rand() % 2);
  }

  return biti;
}


vector<bool> generateNeighbours(int lenght, vector<bool> &biti) {
  // for n dimensions we have n numbers
  // each number has lenght neighbours
  // lenght= number lenght
  vector<bool> neighboursPerDimension;

    for (int j = 0; j < lenght; j++) 
    {
      vector<bool> actNeighbours = biti;
      for(int k=0;k<dimension;k++)
        { if(actNeighbours[lenght*k+j]==false)
            actNeighbours[lenght*k+j]=true;
          else
            actNeighbours[lenght*k+j]=false;
        }
      neighboursPerDimension.insert(neighboursPerDimension.end(),
                                    actNeighbours.begin(), actNeighbours.end());
    }
  
  biti.insert(biti.end(), neighboursPerDimension.begin(),
              neighboursPerDimension.end());

    return biti;
}

// Help for startingIndex: each set has lenght of dimension * requiredBits
vector<double> convertBase(std::vector<bool>& bitstring, int startingIndex = 0)
{
	// Clear the vector first
	while (!converted.empty())
	{
		converted.pop_back();
	}
	for (int i = startingIndex; i < dimension * requiredBits + startingIndex; i = i + requiredBits)
	{
		double convertedNumber = 0;
		for (int j = i; j < i + requiredBits; j++)
		{
			convertedNumber *= 2;
			convertedNumber += bitstring[j];
			// std::cout << candidates[j];
		}
		// std::cout << std::endl;
		convertedNumber = convertedNumber / ((pow(2.000, requiredBits) - 1.000));
		convertedNumber = convertedNumber * (maximum - minimum) + minimum;
		converted.push_back(convertedNumber);
	}

    return converted;
}

long double evaluateFunction(functions& currentFunction)
{
	long double res = 0;
	switch (currentFunction)
	{
	case functions::DeJong:
		res = DeJongFunction();
		break;
	case functions::Schwefel:
		res = SchwefelFunction();
		break;
	case functions::Rastrigin:
		res = RastriginFunction();
		break;
	case functions::Michalewicz:
		res = MichalewiczFunction();
		break;
	default:
		std::cout << "Error at evaluating the function";
		break;
	}
	std::cout << std::fixed << std::setprecision(5) << res;
	return res;
}

vector<bool> firstChoiceHillClimbing(vector<bool> biti)
{
    vector<bool> currentSolution=biti;
    int iterations=0;

    while(iterations<maxIterations)
    { int j=0;
        vector<bool> neighbor= generateNeighbours( requiredBits, currentSolution);
        vector<double> neighborConv=convertBase(neighbor, j);
        vector<double> currentSolConv= convertBase(currentSolution,0);    

        if(neighborConv>currentSolConv)
        { currentSolution=neighbor; j=0;}

		j++;
    cout<< "Iteratie: " << iterations << ", Best: " ;
    for(int k=0;k<currentSolConv.size();k++)
    cout << currentSolConv[k];
    cout<<endl;


     iterations++;   
  }

  return currentSolution;
}

vector<bool> steepestAscentHillClimbing(vector<bool>& biti)
{

    vector<bool> currentSolution = biti; //initial vector
    int iterations = 0;
    
    while (iterations < maxIterations) {
        vector<bool> bestNeighbor = currentSolution;
        vector<double>bestNeighborConv= convertBase(bestNeighbor,0);
        bool foundBetterNeighbor = false;

        for(int i=0; i<requiredBits; ++i)
        {
            vector<bool> neighbor=generateNeighbours(requiredBits,currentSolution);
            vector<double>neighborConv=convertBase(neighbor,i*dimension);
            
        }

        if(!foundBetterNeighbor)
        { break;}

        currentSolution=bestNeighbor;
        cout<<"AAAAAAAAAAAA"<<endl;
        cout << fixed << setprecision(5) << "Iteratie: " << iterations << ", Best: " ;
      for(int k=0;k<bestNeighborConv.size();k++)
        cout << bestNeighborConv[k];
        cout<<endl;
        iterations++;
    }
        return currentSolution;
}



int main()
{
	functions currentFunction;
	currentFunction = functions::Rastrigin;
	setParameters(currentFunction);
    vector<bool>biti;
    generateRandBitNumber(requiredBits,biti);
    vector<bool> FINAL=steepestAscentHillClimbing(biti);

    // converted=convertBase(FINAL,0);
    // cout<< "Solutia finala este: "<<endl;
    // for(int i=0;i<converted.size();i++)
    // cout<<converted[i]<<" ";
    // cout<<endl;

	return 0;
}


void showNeighbours()
{
	int k = 0;
	for (int i = 0; i < bitstring.size(); i++)
	{
		k++;
		std::cout << bitstring[i];
		if (k % requiredBits == 0) std::cout << " ";
		if (k % (requiredBits * dimension) == 0) std::cout << endl;
	}
}

void showConverted()
{
	for (int i = 0; i < converted.size(); i++)
		std::cout << converted[i] << ' ';
	std::cout << endl;
}