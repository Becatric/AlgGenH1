#include <iostream>
#include <fstream> /* output result in a file */
#include <iomanip> /* setprecision */
#include <chrono> /* to calculate the execution time */
#include <time.h> /* for function time(0) used for the seed */
#include <cmath> /* math functions, such as sin, pow, ceil etc*/
#include <vector> /* obviously, vector */
#include <random> /* for mt19937 */
#include <omp.h> /* parallel compiling */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

// Global parameters
int dimension;
double minimum = NULL;
double maximum = NULL;
int precision;
int requiredBits = NULL; // no. of bits required to represent a SINGLE number
ofstream outfile;
char outputLocation[50];
int maxIterations = 20;


std::chrono::time_point<std::chrono::system_clock> start, finish;
std::chrono::duration<double> elapsed_seconds;
mt19937 randomNumber(time(0));

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
	double sum = 0;
	for (int i = 0; i < dimension; i++)
		sum -= sin(converted[i]) * pow(sin((i + 1) * converted[i] * converted[i] / M_PI), 20);
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
		strcpy_s(outputLocation, "DeJong");
		break;
	case functions::Schwefel:
		minimum = -600;
		maximum = 600;
		strcpy_s(outputLocation, "Schwefel");
		break;
	case functions::Rastrigin:
		minimum = -5.12;
		maximum = 5.12;
		strcpy_s(outputLocation, "Rastrigin");
		break;
	case functions::Michalewicz:
		minimum = 0;
		maximum = M_PI; // 3.14
		strcpy_s(outputLocation, "Michalewicz");
		break;
	default:
		std::cout << "Error at setting the function domain";
		break;
	}
	requiredBits = std::ceil(log2((maximum - minimum) * pow(10, precision)));
}

// Generate neighbours - will return a new bool vector, not modifying the source
vector<bool> generateNeighbours(vector<bool> biti) {
	vector<bool> neighboursPerDimension;
	for (int j = 0; j < requiredBits; j++) {
		vector<bool> actNeighbours = biti;
		// Flip the j-st bit of each number
		for (int k = 0; k < dimension; k++) {
			actNeighbours[requiredBits * k + j] = !actNeighbours[requiredBits * k + j];
		}
		// Add the modified neighbors to the vector
		neighboursPerDimension.insert(neighboursPerDimension.end(),
			actNeighbours.begin(), actNeighbours.end());
	}
	return neighboursPerDimension; // Return only the neighbors, not the initial solution
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

//
// Hill Climbing functions
//

// Random bitstring - used only in HC
vector<bool> generateRandBitNumber(int lungime, vector<bool>& biti)
{
	for (int i = 0; i < lungime * dimension; i++)
		biti.push_back(randomNumber() % 2);
	return biti;
}

// First Improvement
double firstChoiceHillClimbing(vector<bool> biti, functions currentFunction)
{
	start = std::chrono::system_clock::now();
	int L = requiredBits * dimension;
	double result = numeric_limits<double>::max();
	for (int iterations = 0; iterations < maxIterations; iterations++)
	{
		vector<bool> currentSolution;
		currentSolution = generateRandBitNumber(requiredBits, currentSolution);
		bool change = true;
		while (change)
		{
			change = false;
			vector<double> currentSolConv = convertBase(currentSolution);
			vector<bool> neighbour = generateNeighbours(currentSolution);
			for (int index = 0; index < requiredBits * L && change == false; index = index + L)
			{
				vector<double> neighbourConv = convertBase(neighbour, index);
				if (evaluateFunction(currentFunction, neighbourConv) < evaluateFunction(currentFunction, currentSolConv))
				{
					change = true;
					vector<bool> aux;
					for (int i = index; i < index + L; i++)
						aux.push_back(neighbour[i]);
					currentSolution = aux;
				}
			}
			if (change)
			{
				currentSolConv = convertBase(currentSolution);
				if (evaluateFunction(currentFunction, currentSolConv) < result)
					result = evaluateFunction(currentFunction, currentSolConv);
			}
		}
	}
	finish = std::chrono::system_clock::now();
	elapsed_seconds = finish - start;
	outfile << "Rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	cout << fixed << setprecision(5) << "Rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	return result;
}

// Best Improvement
double bestChoiceHillClimbing(vector<bool> biti, functions currentFunction)
{
	start = std::chrono::system_clock::now();
	int L = requiredBits * dimension;
	double result = numeric_limits<double>::max();
	for (int iterations = 0; iterations < maxIterations; iterations++)
	{
		vector<bool> currentSolution;
		currentSolution = generateRandBitNumber(requiredBits, currentSolution);
		bool change = true;
		while (change)
		{
			change = false;
			vector<double> currentSolConv = convertBase(currentSolution);
			double doubleCurrentSolution = evaluateFunction(currentFunction, currentSolConv);
			vector<bool> neighbour = generateNeighbours(currentSolution);
			double bestChange = numeric_limits<double>::max();
			int pointer;
			for (int index = 0; index < requiredBits * L; index = index + L)
			{
				vector<double> neighbourConv = convertBase(neighbour, index);
				double doubleNeighbour = evaluateFunction(currentFunction, neighbourConv);
				if (doubleNeighbour < doubleCurrentSolution && doubleNeighbour < bestChange)
				{
					change = true;
					bestChange = doubleNeighbour;
					pointer = index;
				}
			}
			if (change)
			{
				vector<bool> aux;
				for (int i = pointer; i < pointer + L; i++)
					aux.push_back(neighbour[i]);
				currentSolution = aux;
				currentSolConv = convertBase(currentSolution);
				if (evaluateFunction(currentFunction, currentSolConv) < result)
					result = evaluateFunction(currentFunction, currentSolConv);
			}
		}
	}
	finish = std::chrono::system_clock::now();
	elapsed_seconds = finish - start;
	outfile << "Rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	cout << fixed << setprecision(5) << "Rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	return result;
}

double worstChoiceHillClimbing(vector<bool> biti, functions currentFunction)
{
	start = std::chrono::system_clock::now();
	int L = requiredBits * dimension;
	double result = numeric_limits<double>::max();
	for (int iterations = 0; iterations < maxIterations; iterations++)
	{
		vector<bool> currentSolution;
		currentSolution = generateRandBitNumber(requiredBits, currentSolution);
		bool change = true;
		while (change)
		{
			change = false;
			vector<double> currentSolConv = convertBase(currentSolution);
			double doubleCurrentSolution = evaluateFunction(currentFunction, currentSolConv);
			vector<bool> neighbour = generateNeighbours(currentSolution);
			double bestChange = numeric_limits<double>::min();
			int pointer;
			for (int index = 0; index < requiredBits * L; index = index + L)
			{
				vector<double> neighbourConv = convertBase(neighbour, index);
				double doubleNeighbour = evaluateFunction(currentFunction, neighbourConv);
				if (doubleNeighbour < doubleCurrentSolution && doubleNeighbour > bestChange)
				{
					change = true;
					bestChange = doubleNeighbour;
					pointer = index;
				}
			}
			if (change)
			{
				vector<bool> aux;
				for (int i = pointer; i < pointer + L; i++)
					aux.push_back(neighbour[i]);
				currentSolution = aux;
				currentSolConv = convertBase(currentSolution);
				if (evaluateFunction(currentFunction, currentSolConv) < result)
					result = evaluateFunction(currentFunction, currentSolConv);
			}
		}
	}
	finish = std::chrono::system_clock::now();
	elapsed_seconds = finish - start;
	outfile << "Rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	cout << fixed << setprecision(5) << "Rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	return result;
}


// Itearated Simulated Annealing
double simulatedAnnealing(functions& currentFunction, vector<bool>& bitstring)
{
	int l = dimension * requiredBits;
	double candidateResult, neighbourResult, temperature, global = numeric_limits<double>::max();

	for (int i = 0; i < 100; i++)
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

		double t = 1;
		// Start of SA algorithm
		do
		{
			bool no_change = false;
			int contor = 0, contor2 = 0;
			do
			{
				if (!no_change)
				{
					index = 0;
					// Check for best Hamming 1 neighbour
					double bestNeighbour = numeric_limits<double>::max();

					vector<bool>neighbours = generateNeighbours(bitstring);
					bitstring.insert(bitstring.end(), neighbours.begin(), neighbours.end());
					// At index k we have the best neighbour
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
				}

				// Final evaluation of the neighour
				vector<bool> newCandidate;
				for (int k = index; k < index + l; k++)
					newCandidate.push_back(bitstring[k]);

				// If a better Hamming 1 neighbour is not found, try a random one.
				if (index == 0)
				{
					int NumberOfChanges = randomNumber() % (2 * requiredBits - (randomNumber() % 6));
					while (NumberOfChanges)
					{
						int randomIndex = randomNumber() % l;
						newCandidate[randomIndex] = !newCandidate[randomIndex];
						NumberOfChanges--;
					}
				}

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
					contor2++;
				}
			} while (contor < 2 * requiredBits && contor2 < 3);
			t *= 1.7;
			temperature = 1.7 / log2(t + 50); // 70-100
		} while (temperature > 10e-5);

		//cout << fixed << setprecision(5) << "Iteratie: " << i + 1 << ", Best: " << global <<
		//	", local: " << candidateResult << endl;
	}
	//cout << global;
	return global;
}

void compile(functions& currentFunction)
{
	// SIMULATED ANNEALING
	strcat_s(outputLocation, "SA.txt");
	outfile.open(outputLocation, std::ios_base::app);
	outfile << "Results for " << dimension << " dimensions" << '\n' << '\n';
	#pragma omp parallel for
	for (int iteratie = 1; iteratie <= 30; iteratie++)
	{
		start = std::chrono::system_clock::now();
		vector<bool>bitstring;
		double result = simulatedAnnealing(currentFunction, bitstring);
		finish = std::chrono::system_clock::now();
		elapsed_seconds = finish - start;
		outfile << "Iteratie: " << iteratie << ", rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	}
	outfile.close();
}

int main()
{
	functions currentFunction;

	currentFunction = functions::Rastrigin;

	// SA
	/*
	dimension = 5, precision = 7;
	setParameters(currentFunction);
	compile(currentFunction);
	dimension = 10, precision = 3;
	setParameters(currentFunction);
	compile(currentFunction);
	dimension = 30, precision = 1;
	setParameters(currentFunction);
	compile(currentFunction);

	currentFunction = functions::DeJong;
	dimension = 5, precision = 7;
	setParameters(currentFunction);
	compile(currentFunction);
	dimension = 10, precision = 3;
	setParameters(currentFunction);
	compile(currentFunction);
	dimension = 30, precision = 1;
	setParameters(currentFunction);
	compile(currentFunction);

	currentFunction = functions::Schwefel;
	dimension = 5, precision = 7;
	setParameters(currentFunction);
	compile(currentFunction);
	dimension = 10, precision = 3;
	setParameters(currentFunction);
	compile(currentFunction);
	dimension = 30, precision = 1;
	setParameters(currentFunction);
	compile(currentFunction);

	currentFunction = functions::Michalewicz;
	dimension = 5, precision = 7;
	setParameters(currentFunction);
	compile(currentFunction);
	dimension = 10, precision = 3;
	setParameters(currentFunction);
	compile(currentFunction);
	dimension = 30, precision = 1;
	setParameters(currentFunction);
	compile(currentFunction);
	*/

	// HC
	double x;

	// Rastrigin 5D
	dimension = 5, precision = 5;
	setParameters(currentFunction);
	vector<bool> biti;
	biti = generateRandBitNumber(requiredBits, biti);
	
	cout << "First" << endl;
	for (int i = 0; i <= 30; i++)
		x = firstChoiceHillClimbing(biti, currentFunction);
	cout << "Best" << endl;
	for (int i = 0; i <= 30; i++)
		x = bestChoiceHillClimbing(biti, currentFunction);
	cout << "Worst" << endl;
	for (int i = 0; i <= 30; i++)
		x = worstChoiceHillClimbing(biti, currentFunction);
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