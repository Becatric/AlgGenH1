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
int dimension = 10;
double minimum = NULL;
double maximum = NULL;
int precision = 10;
int requiredBits = NULL; // no. of bits required to represent a SINGLE number
std::vector<bool> bitstring; // will have the size of requiredBits * dimension, therefore being able to contain all required numbers for a single function evaluation
std::vector<double> converted; // converted numbers (from the bitstring vector candidates)



std::mt19937 randomNumber(time(0));

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


// add a parameter for the starting position (there will also be the neighbours in bitstring vector
void convertBase(std::vector<bool>& bitstring)
{
	for (int i = 0; i < dimension * requiredBits; i = i + requiredBits)
	{
		long double convertedNumber = 0;
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
	for (int i = 0; i < converted.size(); i++)
		std::cout << converted[i] << std::endl;
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

void simulatedAnnealing(functions& currentFunction)
{
	long double res;
	int t = 0, temperature = 1000000;
	for (int i = 0; i < dimension * requiredBits; i++)
		bitstring.push_back(randomNumber() % 2);
	convertBase(bitstring);
	res = evaluateFunction(currentFunction);
}

int main()
{
	functions currentFunction;
	currentFunction = functions::Rastrigin;
	setParameters(currentFunction);
	simulatedAnnealing(currentFunction);
	return 0;
}