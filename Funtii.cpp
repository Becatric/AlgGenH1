#include<iostream> 
#include<cmath>
#include<vector>
#define PI 3.1415926535897932384626

using namespace std;

double DeJongFuntion(int dimensiuni, vector<double>&x)
{
    double sum=0;

    for (int i = 0; i < dimensiuni; i++)
        sum += x[i] * x[i];

    return sum;
}

double SchwefelFuntion(int dimensiuni, vector<double>&x)
{
    double sum=0;

    for (int i = 0; i < dimensiuni; i++)
        sum -= x[i] * sin(sqrt(abs(x[i])));

    return 418.9829 * dimensiuni + sum;
}

double RastringinFuntion(int dimensiuni, vector<double>&x)
{
    double sum=10*dimensiuni;

    for (int i = 0; i < dimensiuni; i++)
     sum += x[i] * x[i] - 10 * cos(2 * PI * x[i]);

    return sum;
}

double MichalewiczFuntion(int dimensiuni, vector<double>&x)
{ 
    double sum=0;
    int m=10;

    for (int i = 0; i < dimensiuni; i++)
        sum -= sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / PI), 2 * m);

    return sum;
}


int main()
{return 0;}
