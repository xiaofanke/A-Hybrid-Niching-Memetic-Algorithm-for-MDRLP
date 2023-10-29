#pragma once
//#ifndef __MRLP_H_
//#define __MRLP_H_
#define maxNum 9999
#define myRand() rand() % (maxNum + 1) / (double)(maxNum + 1)
#define randAB(a, b) (rand() % (b-a+1))+ a;

#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <string>
#include <map>
#include <time.h>
#include <ilcplex/ilocplex.h>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <future>
#include <chrono>
#include <random>
#include <algorithm>
#include <queue>
#include <set>
#include <stack>

using namespace std;

//typedef IloArray<IloNumVarArray> IloNumVarArray2;
//typedef IloArray<IloBoolVarArray> IloBoolVarArray2;
//typedef IloArray<IloIntVarArray> IloIntVarArray2;
//typedef IloArray<IloArray<IloBoolVarArray>> IloBoolVarArray3;

