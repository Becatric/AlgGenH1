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
int dimension = 2;
double minimum = NULL;
double maximum = NULL;
int precision = 2; // @@@@@@ Cel putin in SA conteaza precizia. Nu stiu daca in HC da, insa poti incerca sa o mai schimbi
// Cu cat precizia e mai mare, cu atat vectorul de biti e mai lung, si, implicit reprezinti mai multe numere
int requiredBits = NULL; // no. of bits required to represent a SINGLE number
int maxIterations = 10;

// WILL BE DELETED
//vector<bool> bitstring; // will have the size of requiredBits * dimension, therefore being able to contain all required numbers for a single function evaluation
//vector<double> converted; // converted numbers (from the bitstring vector candidates)




mt19937 randomNumber(time(0));
// @@@@@@@@@@@@@@@@@@@@ imi genereaza o eroare, nu stiu de ce. Poate pentru ca reqBits e null?
//uniform_int_distribution<int> distribution(1, requiredBits);

enum functions
{
	DeJong = 1,
	Schwefel,
	Rastrigin,
	Michalewicz
};

// Debugging functions - Don't delete them yet, let 'em be
void showNeighbours(vector<bool>& bitstring);
void showConverted(vector<double>& converted);


// The functions themselves
double DeJongFunction(vector<double>& converted)
{
	double sum = 0;
	for (int i = 0; i < dimension; i++)
		sum += converted[i] * converted[i];
	return sum;
}

double SchwefelFunction(vector<double>& converted)
{
	double sum = 0;
	for (int i = 0; i < dimension; i++)
		sum -= converted[i] * sin(sqrt(abs(converted[i])));
	return 418.9829 * dimension + sum;
}

double RastriginFunction(vector<double>& converted)
{
	double sum = 10 * dimension;
	for (int i = 0; i < dimension; i++)
		sum += converted[i] * converted[i] - 10 * cos(2 * M_PI * converted[i]);
	return sum;
}

double MichalewiczFunction(vector<double>& converted)
{
	double sum = 0; int m = 10;
	for (int i = 0; i < dimension; i++)
		sum -= sin(converted[i]) * pow(sin((i + 1) * converted[i] * converted[i] / M_PI), 2 * m);
	return sum;
}



// Based on the current function, set the parameters (minimum, maximum, requiredBits)
void setParameters(functions& currentFunction)
{
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


vector<bool> generateRandBitNumber(int lungime, vector<bool> &biti) 
{
  for (int i = 0; i < lungime; i++)
	  biti.push_back(rand() % 2);
  return biti;
}


// Generate neighbours - they will be added at the final of bitstring, not in a new vector
vector<bool> generateNeighbours(vector<bool> &biti) {
  // for n dimensions we have n numbers
  // each number has <requiredBits> neighbours
  vector<bool> neighboursPerDimension;

    for (int j = 0; j < requiredBits; j++) 
    {
      vector<bool> actNeighbours = biti;
      for(int k=0;k<dimension;k++)
      { 
		  if(actNeighbours[requiredBits*k+j]==false)
             actNeighbours[requiredBits*k+j]=true;
          else
             actNeighbours[requiredBits*k+j]=false;
      }
      neighboursPerDimension.insert(neighboursPerDimension.end(),
                                    actNeighbours.begin(), actNeighbours.end());
    }
	biti.insert(biti.end(), neighboursPerDimension.begin(),
               neighboursPerDimension.end());
    return biti;
}

// Help for startingIndex: each set has lenght of dimension * requiredBits
// Will return a new vector containing the numbers
vector<double> convertBase(std::vector<bool>& bitstring, int startingIndex = 0)
{
	vector<double> converted;
	for (int i = startingIndex; i < dimension * requiredBits + startingIndex; i = i + requiredBits)
	{
		double convertedNumber = 0;
		for (int j = i; j < i + requiredBits; j++)
		{
			convertedNumber *= 2;
			convertedNumber += bitstring[j];
		}
		convertedNumber = convertedNumber / ((pow(2.000, requiredBits) - 1.000));
		convertedNumber = convertedNumber * (maximum - minimum) + minimum;
		converted.push_back(convertedNumber);
	}
    return converted;
}


// Call after convertBase() function to get the effective result
double evaluateFunction(functions& currentFunction, vector<double>& converted)
{
	double res = 0;
	switch (currentFunction)
	{
	case functions::DeJong:
		res = DeJongFunction(converted);
		break;
	case functions::Schwefel:
		res = SchwefelFunction(converted);
		break;
	case functions::Rastrigin:
		res = RastriginFunction(converted);
		break;
	case functions::Michalewicz:
		res = MichalewiczFunction(converted);
		break;
	default:
		std::cout << "Error at evaluating the function";
		break;
	}
	return res;
}

// HillClimbing functions
// @@@@@@@@@@@@@@@@@@@@ Ar putea fi dat parametrul prin referinta sa nu faca o copie?
vector<bool> firstChoiceHillClimbing(vector<bool> biti)
{
	// Nu cred ca are sens sa copiezi intr-un nou vector parametrul. 
	// Nu stiu cum functioneaza functia totusi, poate ai nevoie sa nu
	// schimbi vectorul biti, caz in care e ok.
    vector<bool> currentSolution=biti;
    int iterations=0;

    while(iterations<maxIterations)
    { int j=0; 
	// @@@@@@@@@@@@@@ De ce incepi de la 0? 
	// generateNeighbours da append la vecini. La 0 e candidatul efectiv, 
	// pierzi timp recalculandu-l. Ar trebui sa inceapa la requiredBits*dimension
	// Personal in SA mi-am facut o variabila L = requiredBits*dimension
	// Sa nu calculez de prea multe ori
        vector<bool> neighbor= generateNeighbours(currentSolution);
        vector<double> neighborConv=convertBase(neighbor, j);
        vector<double> currentSolConv= convertBase(currentSolution,0);    

	// @@@@@@@@@@@@ Aici am zis deja ca trebuie un evaluateFunction
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

// @@@@@@@@ idk, prefer sa-mi dai poze cu Hanna decat sa incerc sa inteleg. Prea multa info strica
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
            vector<bool> neighbor=generateNeighbours(currentSolution);
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


// Itearated Simulated Annealing
double simulatedAnnealing(functions& currentFunction, vector<bool>& bitstring)
{
	int l = dimension * requiredBits;
	double candidateResult, neighbourResult, temperature, global = numeric_limits<double>::max();

	for (int i = 0; i < 1000; i++)
	{
		int index = 0;
		temperature = 100;

		// Resetting bitstring
		while (bitstring.size())
			bitstring.pop_back();
		for (int j = 0; j < l; j++)
			bitstring.push_back(randomNumber() % 2);

		vector<double> converted = convertBase(bitstring);
		candidateResult = evaluateFunction(currentFunction, converted);
		if (candidateResult < global)
			global = candidateResult;

		// Start of SA algorithm
		do
		{
			bool no_change = false;
			int contor = 0;
			do
			{
				// Check for best Hamming 1 neighbour
				double bestNeighbour = numeric_limits<double>::max();
				index = 0;
				bitstring = generateNeighbours(bitstring);
				for (int k = l; k <= l * requiredBits; k = k + l)
				{
					converted = convertBase(bitstring, k);
					neighbourResult = evaluateFunction(currentFunction, converted);
					if (neighbourResult < candidateResult && neighbourResult < bestNeighbour)
					{
						index = k;
						bestNeighbour = neighbourResult;
					}
					if (neighbourResult < global)
						global = neighbourResult;
				}
				// If a better neighbour is not found, get a random Hamming 2 one
				if (index == 0)
					index = l * (randomNumber() % requiredBits + 1);

				// Final evaluation of the neighour
				vector<bool> newCandidate;
				for (int k = index; k < index + l; k++)
					newCandidate.push_back(bitstring[k]);
				converted = convertBase(newCandidate);
				neighbourResult = evaluateFunction(currentFunction, converted);

				// First case: a better neighbour is found
				if (neighbourResult < candidateResult)
				{
					//cout << "a";
					bitstring = newCandidate;
					candidateResult = neighbourResult;
					if (candidateResult < global)
						global = candidateResult;
				}
				// Second case: Use a worse neighbour
				else if ((randomNumber() % 1000) / 1000.0 < exp(-abs(neighbourResult - candidateResult) / temperature))
				{
					//cout << "b";
					bitstring = newCandidate;
					candidateResult = neighbourResult;
					if (candidateResult < global)
						global = candidateResult;
					contor++;
				}
				// Third case: Nothing changed, exit do-while and change temperature
				else
				{
					//cout << "c";
					while (bitstring.size() > l)
						bitstring.pop_back();
					no_change = true;
				}
				// Stop when nothing changes or when we reach enough "Second case"
			} while (contor < requiredBits && no_change == false);
			temperature = temperature * 0.9;
		} while (temperature > 10e-8);
		cout << fixed << setprecision(5) << "Iteratie: " << i << ", Best: " << global << ", local: " << candidateResult << endl;
	}
	cout << global;
	return global;
}


int main()
{
	functions currentFunction;
	currentFunction = functions::Rastrigin;
	setParameters(currentFunction);
	
    //vector<bool>biti;
    //generateRandBitNumber(requiredBits,biti);
    //vector<bool> FINAL=steepestAscentHillClimbing(biti);

    // converted=convertBase(FINAL,0);
    // cout<< "Solutia finala este: "<<endl;
    // for(int i=0;i<converted.size();i++)
    // cout<<converted[i]<<" ";
    // cout<<endl;

	vector<bool>bitstring;
	simulatedAnnealing(currentFunction, bitstring);
	return 0;
}


// Debugging functions
void showNeighbours(vector<bool>& bitstring)
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

void showConverted(vector<double>& converted)
{
	for (int i = 0; i < converted.size(); i++)
		std::cout << converted[i] << ' ';
	std::cout << endl;
}