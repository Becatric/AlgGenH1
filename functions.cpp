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
int precision = 5;
int requiredBits = NULL; // no. of bits required to represent a SINGLE number
std::vector<bool> bitstring; // will have the size of requiredBits * dimension, therefore being able to contain all required numbers for a single function evaluation
std::vector<double> converted; // converted numbers (from the bitstring vector)

// Debugging functions - Don't delete them yet, let 'em be
void showNeighbours();
void showConverted();

std::mt19937 randomNumber(time(0));

enum functions
{
	DeJong = 1,
	Schwefel,
	Rastrigin,
	Michalewicz
};

double DeJongFunction()
{
	double sum = 0;

	for (int i = 0; i < dimension; i++)
		sum += converted[i] * converted[i];

	return sum;
}

double SchwefelFunction()
{
	double sum = 0;

	for (int i = 0; i < dimension; i++)
		sum -= converted[i] * sin(sqrt(abs(converted[i])));

	return 418.9829 * dimension + sum;
}

double RastriginFunction()
{
	double sum = 10 * dimension;

	for (int i = 0; i < dimension; i++)
		sum += converted[i] * converted[i] - 10 * cos(2 * M_PI * converted[i]);

	return sum;
}

double MichalewiczFunction()
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

void generateNeighbours(vector<bool>& bitstring) {
	// for n dimensions we have n numbers
	// each number has lenght neighbours
	vector<bool> neighboursPerDimension;

	for (int j = 0; j < requiredBits; j++)
	{
		vector<bool> actNeighbours = bitstring;
		for (int k = 0; k < dimension; k++)
		{
			if (actNeighbours[requiredBits * k + j] == false)
				actNeighbours[requiredBits * k + j] = true;
			else
				actNeighbours[requiredBits * k + j] = false;
		}
		neighboursPerDimension.insert(neighboursPerDimension.end(),
			actNeighbours.begin(), actNeighbours.end());
	}

	bitstring.insert(bitstring.end(), neighboursPerDimension.begin(),
		neighboursPerDimension.end());
}


// Help for startingIndex: each set has lenght of dimension * requiredBits
void convertBase(std::vector<bool>& bitstring, int startingIndex = 0)
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
}

double evaluateFunction(functions& currentFunction)
{
	double res = 0;
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
	return res;
}

double simulatedAnnealing(functions& currentFunction)
{
	int l = dimension * requiredBits;
	double candidateResult, neighbourResult, temperature, global=9999999;
	for (int i = 0; i < 1000; i++)
	{
		int index = 0;
		temperature = 100;
		while (bitstring.size())
			bitstring.pop_back();
		for (int j = 0; j < l; j++)
			bitstring.push_back(randomNumber() % 2);
		convertBase(bitstring);
		candidateResult = evaluateFunction(currentFunction);
		if (candidateResult < global)
			global = candidateResult;
		do
		{
			bool no_change = false;
			int contor = 0;
			do
			{
				int bestNeighbour = 999999999;
				index = 0;
				generateNeighbours(bitstring);
				for (int k = l; k <= l * requiredBits; k = k + l)
				{
					convertBase(bitstring, k);
					neighbourResult = evaluateFunction(currentFunction);
					if (neighbourResult < candidateResult && neighbourResult < bestNeighbour)
					{
						index = k;
						bestNeighbour = neighbourResult;
					}
					if (neighbourResult < global)
						global = neighbourResult;
				}
				if (index == 0)
					index = l * (randomNumber() % requiredBits + 1);

				vector<bool> newCandidate;
				for (int k = index; k < index + l; k++)
					newCandidate.push_back(bitstring[k]);
				convertBase(newCandidate);
				neighbourResult = evaluateFunction(currentFunction);

				if (neighbourResult < candidateResult)
				{
					cout << "a";
					bitstring = newCandidate;
					candidateResult = neighbourResult;
					if (candidateResult < global)
						global = candidateResult;
				}
				else if ((randomNumber() % 1000) / 1000.0 < exp(-abs(neighbourResult - candidateResult) / temperature))
				{
					cout << "b";
					bitstring = newCandidate;
					candidateResult = neighbourResult;
					if (candidateResult < global)
						global = candidateResult;
					contor++;
				}
				else
				{
					cout << "c";
					while (bitstring.size() > l)
						bitstring.pop_back();
					no_change = true;
				}
			// Stop when nothing changes or when we reach enough "else if"
			} while (contor < requiredBits && no_change == false);
			temperature = temperature * 0.9;
		} while (temperature > 10e-8);
		std::cout << fixed << setprecision(5) << "Iteratie: " << i << ", Best: " << global << ", local: " << candidateResult << endl;
	}
	std::cout << global;
	return global;
}

int main()
{
	functions currentFunction;
	currentFunction = functions::Rastrigin;
	setParameters(currentFunction);
	simulatedAnnealing(currentFunction);
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