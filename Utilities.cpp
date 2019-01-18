#include "Utilities.hpp"

using namespace std;

string itostr(int number)
{
  stringstream ss;//create a stringstream
  ss << number;//add number to the stream
  return ss.str();//return a string with the contents of the stream
}

// double getRandom()// get a die number \in [0, 1]
// {
//   return double (rand())/double(RAND_MAX);
// }

