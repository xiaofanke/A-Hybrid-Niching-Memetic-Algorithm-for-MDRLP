#pragma once
#include "head.h"
#include "problem.h"

class IndividualClass {
public:
	vector<int> sequence; // ���н���

	vector<vector<int>> layout; // ���н�
	vector<double> location; // ����ֵ
	vector<int> indexR; // �豸���ڵ���
	double offset;

	double obj_offset;
	double obj_final;

	int node_index[51]; // ��¼�豸�������е�λ��

	// ��ʼ��һ���⣬����
	void init(ProblemClass problem);

	// �����˵�
	void adjustBreakpoint(ProblemClass problem);

	// ����
	double getObj_offset(int& objNum, ProblemClass problem, int num);


	double getObj_final(ProblemClass problem, int num);


	//����Ŀ�꺯��ֵ
	void write(ProblemClass problem, int num);
};