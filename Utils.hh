#ifndef Utils_hh
#define Utils_hh

using namespace std;

/////////////////// random numbers //////////////////////

// Returns a random number between l and u.
int randomize(int l, int u);


//////////////////// auxiliary types for identifiers ////////////////////

typedef int   Player_Id;
typedef int Starship_Id;
typedef int  Missile_Id;


/////////////////// positions, directions  //////////////

#include <complex>

// If needed, consult the complex<T> template class from the STL.
typedef complex<int> Pos;
typedef complex<int> Dir;


/////////////////////////////////////
//                                 //
//    ALL ALLOWED DIRECTIONS       //
//                                 //
constexpr Dir UP        = {-1, 1}; //
constexpr Dir DEFAULT   = { 0, 1}; //
constexpr Dir DOWN      = { 1, 1}; //
//                                 //
constexpr Dir SLOW_UP   = {-1, 0}; //
constexpr Dir SLOW      = { 0, 0}; //
constexpr Dir SLOW_DOWN = { 1, 0}; //
//                                 //
constexpr Dir FAST_UP   = {-1, 2}; //
constexpr Dir FAST      = { 0, 2}; //
constexpr Dir FAST_DOWN = { 1, 2}; //
//                                 //
/////////////////////////////////////

/**
 * Returns the first component (rows)
 */
int  first(complex<int> x);


/**
 * Returns the second component (columns)
 */
int second(complex<int> x);


/****************************************************/
/***************** STOP READING HERE ****************/
/****************************************************/


inline bool dir_ok(Dir d) {
  return
    d == UP        or
    d == DEFAULT   or
    d == DOWN      or
    d == SLOW_UP   or
    d == SLOW      or
    d == SLOW_DOWN or
    d == FAST_UP   or
    d == FAST      or
    d == FAST_DOWN;
}

inline int  first(complex<int> x) { return x.real(); }
inline int second(complex<int> x) { return x.imag(); }

#include <map>
#include <set>
#include <vector>

struct Cmp {
  bool operator()(Dir a, Dir b) const {
    return make_pair(first(a), second(a)) < make_pair(first(b), second(b));
  }
};

extern const map<Dir, vector<Dir>, Cmp> dir2all;
  
/////////////////// random numbers //////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <queue>
#include <limits>
#include <cassert>
#include <cstdlib>
#include <cerrno>
#include <time.h>
#include <sys/time.h>

class Random_generator {

  static constexpr int RANDOMIZE_MODULUS    = 2147483647;
  static constexpr int RANDOMIZE_MULTIPLIER = 48271;
  static constexpr int RANDOMIZE_DEFAULT    = 123456789;  // Initial seed.
  static constexpr int RANDOMIZE_MAX        = RANDOMIZE_MODULUS - 2;
      
  long rnd_seed;

  // Sets random seed.
  void srandomize(int s = RANDOMIZE_DEFAULT) {
    if (s < 0) s = -s;
    rnd_seed = 1 + s % (RANDOMIZE_MODULUS - 1);
  }

  friend int main (int argc, char** argv);

  
public:

  // Returns random number between 0 and RANDOMIZE_MAX.
  int randomize(void) {
    constexpr long Q = RANDOMIZE_MODULUS / RANDOMIZE_MULTIPLIER;
    constexpr long R = RANDOMIZE_MODULUS % RANDOMIZE_MULTIPLIER;
    long t;

    t = RANDOMIZE_MULTIPLIER * (rnd_seed % Q) - R * (rnd_seed / Q);
    if (t > 0) 
      rnd_seed = t;
    else 
      rnd_seed = t + RANDOMIZE_MODULUS;
    return rnd_seed - 1;
  }
  
};


extern Random_generator random_generator;

inline int randomize() {
  return random_generator.randomize();
}
  
// Returns a random number between l and u.
inline int randomize(int l, int u) {
  return l + randomize() % (u-l+1);

}

// Returns a random permutation of {0, 1, ..., n-1}.
inline vector<int> random_permutation(int n) {
  vector<int> v(n);
  for (int i = 0; i < n; ++i) v[i] = i;
  for (int i = 0; i < n; ++i) swap(v[i], v[randomize(i, n-1)]);
  return v;
}

////////////////////////////////////////////////////////////////

#define contains count

////////////////////////////////////////////////////////////////

#define my_assert(X)                                                              \
  if (not (X)) {							\
       cerr << "Failed assertion '" << #X << "' in file "			  \
                                    << __FILE__ << ", line " << __LINE__ << "\n"; \
       exit(1);                                                                   \
    }

#endif // Utils_hh
