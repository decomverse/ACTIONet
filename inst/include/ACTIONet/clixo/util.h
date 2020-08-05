#ifndef UTIL
#define UTIL

#include <string>
#include <vector>

using namespace std;

class Utils {
 public:
  template <typename T>
    static void printVector(const vector<T> & vec) {
    if (vec.size() != 0) {
      typename vector<T>::const_iterator it = vec.begin();
      cout << *it;
      ++it;
      for (; it != vec.end(); ++it) {
	cout << " " << *it;
      }
    }
    cout << endl;
    return;
  }

  inline static void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ") {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);
    
    while (string::npos != pos || string::npos != lastPos) {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
  }

  template <typename T>
    inline static bool insertInOrder(vector<T> & sortedVector, const T & newElement) {
    if ((sortedVector.size() == 0) || (newElement > sortedVector.back())) {
      sortedVector.push_back(newElement);
      return true;
    }
    if (newElement < sortedVector[0]) {
      sortedVector.insert(sortedVector.begin(), newElement);
      return true;
    }
    unsigned minPos = 0;
    unsigned maxPos = sortedVector.size() - 1;
    while ((maxPos - minPos) > 1) {
      unsigned testPos = (minPos + maxPos) / 2;
      T testVal = sortedVector[testPos];
      if (newElement > testVal) {
	minPos = testPos;
      } else if (newElement < testVal) {
	maxPos = testPos;
      } else if (newElement == testVal) {
	return false;
      }
    }
    
    if ((newElement == sortedVector[minPos]) || (newElement == sortedVector[maxPos])) {
      return false;
    }
    //cout << "Inserting " << newElement << "\t" << sortedVector[minPos] << "\t" << sortedVector[maxPos] << endl;
    typename vector<T>::iterator it = sortedVector.begin();
    it += maxPos;
    //printVector(sortedVector);
    sortedVector.insert(it, newElement);
    //printVector(sortedVector);
    return true;
  }

  template <typename T>
    inline static bool eraseInOrder(vector<T> & sortedVector, const T & toErase) {

    signed long minPos = 0;
    signed long maxPos = sortedVector.size() - 1;
    while (maxPos >= minPos) {
      signed long testPos = (minPos + maxPos) / 2;
      T testVal = sortedVector[testPos];
      if (toErase > testVal) {
	minPos = testPos + 1;
      } else if (toErase < testVal) {
	maxPos = testPos - 1;
      } else {
	typename vector<T>::iterator it = sortedVector.begin();
	it += testPos;
	sortedVector.erase(it);
	return true;
      }
    }
    return false;
  }

  template <typename T>
    inline static bool elementExists(vector<T> & sortedVector, const T & toCheck) {
    
    signed long minPos = 0;
    signed long maxPos = sortedVector.size() - 1;
    while (maxPos >= minPos) {
      signed long testPos = (minPos + maxPos) / 2;
      T testVal = sortedVector[testPos];
      if (toCheck > testVal) {
	minPos = testPos + 1;
      } else if (toCheck < testVal) {
	maxPos = testPos - 1;
      } else {
	return true;
      }
    }
    return false;
  }

};

#endif // UTIL
