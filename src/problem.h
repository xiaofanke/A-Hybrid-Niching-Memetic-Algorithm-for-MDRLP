#pragma once
#include "head.h"

class ProblemClass
{
public:

	int facility_num; // 设备个数	
	int row_num; // 设备个数

	vector<double> widths; //设备间最小间隙
	double Cs; //行间距
	vector<vector<double>> minSpaces; //设备间最小间隙
	vector<vector<double>> materials; //设备间最小间隙	

	ProblemClass();
	ProblemClass(char* file_dir);

};
