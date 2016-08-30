// Shape.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <iostream>
#include <math.h>    
  

using namespace std;

long int f(long int N)
{

  long int a=0, b=1 , c ;

  if( N < 1) { 
    cout << "Enter positive number and N should be greater than zero " << endl; 
    return 0;
  }

  for (int i = 2; i <= N; i++)
    {
      c = a + b;
      a = b;
      b = c;
    }
  return a;
}

int main(int argc, char* argv[])
{
  long int N; 
  long int ans ; 
  cout << "Enter the value " << endl; 

  cin >> N; 

  ans =  f(N); 

  cout << "The Nth value is  " << ans << endl; 
  return 0;
}

