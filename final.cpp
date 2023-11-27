#include <iostream>
#include <fstream> /* output result in a file */
#include <iomanip> /* setprecision */
#include <chrono> /* to calculate the execution time */
#include <time.h> /* for function time(0) used for the seed */
#include <cmath> /* math functions, such as sin, pow, ceil etc*/
#include <vector> /* obviously, vector */
#include <random> /* for mt19937 */
#include <omp.h> /* parallel */
#include <cstring> /* VSC */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

// Global parameters
int dimension;
double minimum;
double maximum;
int precision;
int requiredBits; // no. of bits required to represent a SINGLE number
ofstream resultFile, timeFile;
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
void showBitstring(vector<bool>& bitstring);
void showConverted(vector<double>& converted);


// The functions themselves
double DeJongFunction(const vector<double>& converted)
{
	double sum = 0;
	for (int i = 0; i < dimension; i++)
		sum += converted[i] * converted[i];
	return sum;
}

double SchwefelFunction(const vector<double>& converted)
{
	double sum = 0;
	for (int i = 0; i < dimension; i++)
		sum -= converted[i] * sin(sqrt(abs(converted[i])));
	return 418.9829 * dimension + sum;
}

double RastriginFunction(const vector<double>& converted)
{
	double sum = 10 * dimension;
	for (int i = 0; i < dimension; i++)
		sum += converted[i] * converted[i] - 10 * cos(2 * M_PI * converted[i]);
	return sum;
}

double MichalewiczFunction(const vector<double>& converted)
{
	double sum = 0;
	for (int i = 0; i < dimension; i++)
		sum -= sin(converted[i]) * pow(sin((i + 1) * converted[i] * converted[i] / M_PI), 20);
	return sum;
}



// Based on the current function, set the parameters (minimum, maximum, requiredBits)
void setParameters(const functions& currentFunction)
{
	switch (currentFunction)
	{
	case functions::DeJong:
		minimum = -5.12;
		maximum = 5.12;
		resultFile.open("DeJongResult.txt", std::ios_base::app); // append mode
		timeFile.open("DeJongTime.txt", std::ios_base::app);
		break;
	case functions::Schwefel:
		minimum = -600;
		maximum = 600;
		resultFile.open("SchwefelResult.txt", std::ios_base::app);
		timeFile.open("SchwefelTime.txt", std::ios_base::app);
		break;
	case functions::Rastrigin:
		minimum = -5.12;
		maximum = 5.12;
		resultFile.open("RastriginResult.txt", std::ios_base::app);
		timeFile.open("RastriginTime.txt", std::ios_base::app);
		break;
	case functions::Michalewicz:
		minimum = 0;
		maximum = M_PI;
		resultFile.open("MichalewiczResult.txt", std::ios_base::app);
		timeFile.open("MichalewiczTime.txt", std::ios_base::app);
		break;
	default:
		std::cout << "Error at setting the function domain";
		break;
	}
	requiredBits = std::ceil(log2((maximum - minimum) * pow(10, precision)));
}


void writeToFile(double result, double time)
{
	resultFile << result << '\n';
	timeFile << elapsed_seconds.count() << '\n';
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
double evaluateFunction(const functions& currentFunction, const vector<double>& converted)
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

// First Improvement
double firstChoiceHillClimbing(const functions& currentFunction)
{
	start = std::chrono::system_clock::now();
	int L = requiredBits * dimension;
	double result = numeric_limits<double>::max();
	for (int iterations = 0; iterations < maxIterations; iterations++)
	{
		vector<bool> currentSolution;
		for (int j = 0; j < L; j++)
			currentSolution.push_back(randomNumber() % 2);
		bool change = true;
		while (change)
		{
			change = false;
			vector<double> currentSolConv = convertBase(currentSolution);
			//vector<bool> neighbour = generateNeighbours(currentSolution);
			for (int index = 0; index < L && change == false; index = index + 1)
			{
				currentSolution[index] = !currentSolution[index];
				vector<double> neighbourConv = convertBase(currentSolution);
				if (evaluateFunction(currentFunction, neighbourConv) < evaluateFunction(currentFunction, currentSolConv))
					change = true;
				else
					currentSolution[index] = !currentSolution[index];
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
	writeToFile(result, elapsed_seconds.count());
	cout << fixed << setprecision(precision) << "Rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	return result;
}

// Best Improvement
double bestChoiceHillClimbing(const functions& currentFunction)
{
	start = std::chrono::system_clock::now();
	int L = requiredBits * dimension;
	double result = numeric_limits<double>::max();
	for (int iterations = 0; iterations < maxIterations; iterations++)
	{
		vector<bool> currentSolution;
		for (int j = 0; j < L; j++)
			currentSolution.push_back(randomNumber() % 2);
		bool change = true;
		while (change)
		{
			change = false;
			vector<double> currentSolConv = convertBase(currentSolution);
			double doubleCurrentSolution = evaluateFunction(currentFunction, currentSolConv);
			double bestChange = numeric_limits<double>::max();
			int pointer;
			for (int index = 0; index < L; index++)
			{
				currentSolution[index] = !currentSolution[index];
				vector<double> neighbourConv = convertBase(currentSolution);
				double doubleNeighbour = evaluateFunction(currentFunction, neighbourConv);
				if (doubleNeighbour < doubleCurrentSolution && doubleNeighbour < bestChange)
				{
					change = true;
					bestChange = doubleNeighbour;
					pointer = index;
				}
				currentSolution[index] = !currentSolution[index];
			}
			if (change)
			{
				currentSolution[pointer] = !currentSolution[pointer];
				currentSolConv = convertBase(currentSolution);
				if (evaluateFunction(currentFunction, currentSolConv) < result)
					result = evaluateFunction(currentFunction, currentSolConv);
			}
		}
	}
	finish = std::chrono::system_clock::now();
	elapsed_seconds = finish - start;
	writeToFile(result, elapsed_seconds.count());
	cout << fixed << setprecision(precision) << "Rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	return result;
}

double worstChoiceHillClimbing(const functions& currentFunction)
{
	start = std::chrono::system_clock::now();
	int L = requiredBits * dimension;
	double result = numeric_limits<double>::max();
	for (int iterations = 0; iterations < maxIterations; iterations++)
	{
		vector<bool> currentSolution;
		for (int j = 0; j < L; j++)
			currentSolution.push_back(randomNumber() % 2);
		bool change = true;
		while (change)
		{
			change = false;
			vector<double> currentSolConv = convertBase(currentSolution);
			double doubleCurrentSolution = evaluateFunction(currentFunction, currentSolConv);
			double bestChange = -9999999;
			int pointer;
			for (int index = 0; index < L; index++)
			{
				currentSolution[index] = !currentSolution[index];
				vector<double> neighbourConv = convertBase(currentSolution);
				double doubleNeighbour = evaluateFunction(currentFunction, neighbourConv);
				if (doubleNeighbour < doubleCurrentSolution && doubleNeighbour > bestChange)
				{
					change = true;
					bestChange = doubleNeighbour;
					pointer = index;
				}
				currentSolution[index] = !currentSolution[index];
			}
			if (change)
			{
				currentSolution[pointer] = !currentSolution[pointer];
				currentSolConv = convertBase(currentSolution);
				if (evaluateFunction(currentFunction, currentSolConv) < result)
					result = evaluateFunction(currentFunction, currentSolConv);
			}
		}
	}
	finish = std::chrono::system_clock::now();
	elapsed_seconds = finish - start;
	writeToFile(result, elapsed_seconds.count());
	cout << fixed << setprecision(precision) << "Rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	return result;
}


// Itearated Simulated Annealing
double simulatedAnnealing(const functions& currentFunction, vector<bool>& bitstring)
{
	int l = dimension * requiredBits;
	double candidateResult, neighbourResult, temperature, global = numeric_limits<double>::max();

	for (int i = 0; i < 100; i++)
	{
		bool foundBetter;
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
				foundBetter = false;
				if (!no_change)
				{
					index = 0;
					double bestNeighbour = numeric_limits<double>::max();

					// Check for best Hamming 1 neighbour
					for (int k = 0; k < l; k++)
					{
						bitstring[k]=!bitstring[k];
						converted = convertBase(bitstring);
						neighbourResult = evaluateFunction(currentFunction, converted);
						if (neighbourResult < candidateResult && neighbourResult < bestNeighbour)
						{
							foundBetter = true;
							index = k;
							bestNeighbour = neighbourResult;
						}
						if (neighbourResult < global)
							global = neighbourResult;
						bitstring[k]=!bitstring[k];
					}
				}
				vector<bool> newCandidate = bitstring;
				if(foundBetter)
					bitstring[index]=!bitstring[index];
				// If a better Hamming 1 neighbour is not found, try a random one.
				else
				{
					int NumberOfChanges = randomNumber() % 6;
					while (NumberOfChanges)
					{
						int randomIndex = randomNumber() % l;
						newCandidate[randomIndex] = !newCandidate[randomIndex];
						NumberOfChanges--;
					}
				}

				converted = convertBase(bitstring);
				neighbourResult = evaluateFunction(currentFunction, converted);

				// First case: a better neighbour is found
				if (neighbourResult < candidateResult)
				{
					//cout << "a";
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
					no_change = true;
					contor2++;
				}
			} while (contor < 2 * requiredBits && contor2 < 3);
			t *= 1.7;
			temperature = 1.7 / log2(t); // 70-100
		} while (temperature > 10e-5);

		//cout << fixed << setprecision(precision) << "Iteratie: " << i + 1 << ", Best: " << global <<
		//	", local: " << candidateResult << endl;
	}
	//cout << global;
	return global;
}

void compile(const functions& currentFunction)
{
	resultFile << "Simulated Annealing results for " << dimension << " dimensions" << '\n' << '\n';
	timeFile << "Simulated Annealing results for " << dimension << " dimensions" << '\n' << '\n';
	#pragma omp parallel for
	for (int iteratie = 1; iteratie <= 30; iteratie++)
	{
		start = std::chrono::system_clock::now();
		vector<bool>bitstring;
		double result = simulatedAnnealing(currentFunction, bitstring);
		finish = std::chrono::system_clock::now();
		elapsed_seconds = finish - start;
		writeToFile(result, elapsed_seconds.count());
		cout << fixed << setprecision(precision) << "Iteratie: " << iteratie << ", rezultat: " << result << ", durata executie: " << elapsed_seconds.count() << '\n';
	}
}

int main()
{
	functions currentFunction;
	currentFunction = functions::Rastrigin;
	// 5D, precision 10 // 10D precision 7 // 30D precision 2
	dimension = 30, precision = 2;
	setParameters(currentFunction);
	
	
	double x;
	resultFile << '\n' << "First Improvement results for " << dimension << " dimensions" << '\n' << '\n';
	timeFile << '\n' << "First Improvement results for " << dimension << " dimensions" << '\n' << '\n';
	for (int i = 0; i <= 30; i++)
		x = firstChoiceHillClimbing(currentFunction);
	resultFile << '\n' << "Best Improvement results for " << dimension << " dimensions" << '\n' << '\n';
	timeFile << '\n' << "Best Improvement results for " << dimension << " dimensions" << '\n' << '\n';
	for (int i = 0; i <= 30; i++)
		x = bestChoiceHillClimbing(currentFunction);
	resultFile << '\n' << "Worst Improvement results for " << dimension << " dimensions" << '\n' << '\n';
	timeFile << '\n' << "Worst Improvement results for " << dimension << " dimensions" << '\n' << '\n';
	for (int i = 0; i <= 30; i++)
		x = worstChoiceHillClimbing(currentFunction);
	
	// Simulated Annealing
	compile(currentFunction);

	resultFile.close();
	timeFile.close();
}


// Debugging functions
void showBitstring(const vector<bool>& bitstring)
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

void showConverted(const vector<double>& converted)
{
	for (int i = 0; i < converted.size(); i++)
		std::cout << converted[i] << ' ';
	std::cout << endl;
}