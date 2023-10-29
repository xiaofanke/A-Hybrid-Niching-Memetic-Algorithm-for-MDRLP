#pragma once
#include "head.h"
#include "problem.h"

class IndividualClass {
public:
	vector<int> sequence; // 序列解码

	vector<vector<int>> layout; // 序列解
	vector<double> location; // 坐标值
	vector<int> indexR; // 设备所在的行
	double offset;

	double obj_offset;
	double obj_final;

	int node_index[51]; // 记录设备在序列中的位置

	// 初始化一个解，编码
	void init(ProblemClass problem);

	// 调整端点
	void adjustBreakpoint(ProblemClass problem);

	// 解码
	double getObj_offset(int& objNum, ProblemClass problem, int num);


	double getObj_final(ProblemClass problem, int num);


	//计算目标函数值
	void write(ProblemClass problem, int num);
};