/******************************************************************************

                              Online C++ Compiler.
               Code, Compile, Run and Debug C++ program online.
Write your code in this editor and press "Run" button to compile and execute it.

*******************************************************************************/


using namespace std;

#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <algorithm>


//To Do: MongoDB Routing


//mean of vector values
double VecMean (vector < double > x) {

  int n = x.size ();
  double sum = 0.0;
  
  for (int i = 0; i <= n - 1; i++)
    sum += x[i];

  double xbar = sum / n;
  return xbar;
}


//standard normal N(0,1) density
double f (double x) {

  double pi = 4.0 * atan (1.0);
  return exp (-x * x * 0.5) / sqrt (2 * pi);
}


//Boole's Rule
double Boole (double StartPoint, double EndPoint, int n) {
    
  vector < double >X (n + 1, 0.0);
  vector < double >Y (n + 1, 0.0);
  double delta_x = (EndPoint - StartPoint) / double (n);

  for (int i = 0; i <= n; i++) {
      X[i] = StartPoint + i * delta_x;
      Y[i] = f (X[i]);
    }

  double sum = 0;

  for (int t = 0; t <= (n - 1) / 4; t++) {
      int ind = 4 * t;
      sum += (1 / 45.0) * (14 * Y[ind] + 64 * Y[ind + 1] + 24 * Y[ind + 2] + 64 * Y[ind + 3] + 14 * Y[ind + 4]) * delta_x;

    }
    
  return sum;
}


//standard normal N(0, 1) cumulative distribution function(CDF) by Boole's Rule
double N (double x) {
    
  return Boole (-10.0, x, 240);
}


//Black Scholes closed form price calculations
double BlackScholes (double S, double K, double V, double T, double R, double Q, char PutCall) {

  double d1 = (log (S / K) + (R - Q + 0.5 * V * V) * T) / V / sqrt (T);
  double d2 = d1 - V * sqrt (T);
  double Call = S * exp (-Q * T) * N (d1) - K * exp (-R * T) * N (d2);

  if (PutCall == 'C')
    return Call;

  else
    return Call - S * exp (-Q * T) + K * exp (-R * T);
}


int main () {
  
  //seed setter for random (To Do: accommodate FMP API Options Data here instead of random)
  srand (time (0));  
  //spot price setter
  double S = 100.0;  
  //strike price setter
  double K = 100.0;  
  //years maturity setter
  double T = 1;     
  //interest rate setter
  double R = 0.05;  
  //dividend yields setter
  double Q = 0;     
  //volatility rate setter
  double V = 0.2;   
  //number of simulations setter (for 1 million)
  int Nsims = 1e6; 
  //random standard normal setter
  double Z;   

  //terminal price S(T) initialize
  vector < double > ST (Nsims, 0.0); 
  //call payoff initialize
  vector < double > ST_K (Nsims, 0.0); 
  //put payoff setter
  vector < double > K_ST (Nsims, 0.0);

  double u1, u2;
  double pi = 3.141592653589793;

  for (int i = 0; i <= Nsims - 1; i++) {
      //independent uniform random variables
      u1 = ((double) rand () / ((double) (RAND_MAX) + (double) (1)));
      u2 = ((double) rand () / ((double) (RAND_MAX) + (double) (1)));
      //floor u1 to avoid errors with log function
      u1 = max (u1, 1.0e-10);
      //Z ~ N(0,1) by Box-Muller transformation
      Z = sqrt (-2.0 * log (u1)) * sin (2 * pi * u2);
      //terminal price simulated S(T)
      ST[i] = S * exp ((R - Q - 0.5 * V * V) * T + V * sqrt (T) * Z); 
      //call payoff
      ST_K[i] = max (ST[i] - K, 0.0); 
      //put payoff
      K_ST[i] = max (K - ST[i], 0.0); 
    }

  //simulated call & put prices as discounted average of TP
  double BSCallSim = exp (-R * T) * VecMean (ST_K);
  double BSPutSim = exp (-R * T) * VecMean (K_ST);

  //closed form prices
  double BSCall = BlackScholes (S, K, V, T, R, Q, 'C');
  double BSPut = BlackScholes (S, K, V, T, R, R, 'P');

  //margin or errors
  double CallError = BSCall - BSCallSim;
  double PutError = BSPut - BSPutSim;

  //output format
  cout << setprecision (4) << fixed;
  cout << "From " << Nsims << " simulations: " << endl;
  cout << " " << endl;
  cout << "Michael Kim's Monte Carlo Simulation Results are as Following " << endl;
  cout << "Column 1: Call Price ($), Column 2: Put Price ($) " << endl;
  cout << "----------------------------------" << endl;
  cout << "Mean Simulated Prices: " << BSCallSim << " " << BSPutSim << endl;
  cout << "Mean Closed Form Prices: " << BSCall << " " << BSPut << endl;
  cout << "Mean Margin of Error: " << CallError << " " << PutError << endl;
  cout << "----------------------------------" << endl;
  cout << " " << endl;

  system ("PAUSE");
}