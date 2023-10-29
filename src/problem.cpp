#pragma once
#include "problem.h"

ProblemClass::ProblemClass(char* file_dir)
{
	int tmpi;
	double tmpd;

	fstream Stream(file_dir, ios::in);
	char str[200];

	//设备个数	
	Stream >> tmpi;
	facility_num = tmpi;

	//布局空间行数	
	Stream >> tmpi;
	row_num = tmpi;

	//布局设备间间隙	
	Stream >> tmpd;
	Cs = tmpd;

	double c = 0;
	widths.resize(facility_num + 1);

	for (int j = 1; j < facility_num + 1; j++) {
		Stream >> tmpd;
		widths[j] = tmpd;
	}

	//设备间间隙
	minSpaces.resize(facility_num + 1);
	for (int j = 1; j < facility_num + 1; j++) {
		minSpaces[j].resize(facility_num + 1);
		for (int k = 1; k < facility_num + 1; k++) {
			Stream >> tmpd;
			minSpaces[j][k] = tmpd;
		}
	}

	//设备间物料流
	materials.resize(facility_num + 1);
	for (int j = 1; j < facility_num + 1; j++) {
		materials[j].resize(facility_num + 1);
		for (int k = 1; k < facility_num + 1; k++) {
			Stream >> tmpd;
			materials[j][k] = tmpd;
		}
	}

	Stream.close();
}