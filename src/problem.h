#pragma once
#include "head.h"

class ProblemClass
{
public:

	int facility_num; // �豸����	
	int row_num; // �豸����

	vector<double> widths; //�豸����С��϶
	double Cs; //�м��
	vector<vector<double>> minSpaces; //�豸����С��϶
	vector<vector<double>> materials; //�豸����С��϶	

	ProblemClass();
	ProblemClass(char* file_dir);

};
