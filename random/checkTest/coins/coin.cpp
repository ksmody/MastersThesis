// Shape.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <iostream>
#include <math.h>    
  

using namespace std;

void Coins(long int N, long int& nNumPenny, long int& nNumNickel, long int& nNumQuarter)
{ 
  if (N < 0) { 
    cout << "Enter positive number " << endl; 
    return;
  }
  long int amount;  

  nNumQuarter = N/25 ; 

  amount = N - (nNumQuarter*25); 

  nNumNickel = amount/ 5; 

  amount = amount - (nNumNickel*5);

  nNumPenny = amount ; 

  cout << "The least number of coins are " << nNumQuarter + nNumNickel + nNumPenny << " " << "which consists of " << endl;

  cout << "The number of Quater coins " << nNumQuarter << endl;
  cout << "The number of Nickel coins  " << nNumNickel << endl;
  cout << "The number of pennies  " << nNumPenny << endl; 

}

int main(int argc, char* argv[])
{
  long int N, Penny = 0, Nickel=0, Quarter = 0; 

  cout << "Enter the value in cents " << endl; 

  cin >> N; 

  Coins(N, Penny, Nickel, Quarter) ; 

	
  return 0;
}

