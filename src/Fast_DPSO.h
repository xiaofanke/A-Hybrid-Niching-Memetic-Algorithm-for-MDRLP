#pragma once
#include "head.h"
#include "Individual.h"
#define MAXDIM   60	// max dimension
#define MAXSIZE  200    // max swarm size
#define TBSIZE 15 // Tabu list size

class Fast_DPSOClass {
public:
	int popNum; //��Ⱥ��С
	int maxTime; //������ʱ��
	double WR; //�������
	double CR1; //�������
	double CR2; //�������
	int num_eval; //��Ӧ�����۴���

	vector<IndividualClass> pops; // ������Ⱥ
	vector<IndividualClass> allBest; // �������Ž�

	vector<IndividualClass> pBest; // ÿ�����ӵ���ʷ���Ž�
	vector<IndividualClass> gBest; // ÿ��С�������������Ž�


	vector<IndividualClass> bestSolution;
	vector<IndividualClass> allSolution;


	queue<int> tabu_list; // ���ɱ�
	int tabu_index[MAXSIZE] = { 0 }; // ��¼��ǰ�����Ƿ񱻽��ɣ���߲�ѯ�ٶ�
	
	int nsize; // �ھӵĸ���
	int num_projs; // nh ��ϣ��������
	int neighbor[MAXSIZE][MAXSIZE];

	int ls_swap[1500][2]; // 2top�ֲ��������п�ѡ�������豸����λ�á�
	vector<int> ls_index; // ���п�ѡ�������豸����λ�õ�����

	Fast_DPSOClass(int _popNum, int _maxTime, double _wr, double _cr1, double _cr2, int _num_projs) {
		popNum = _popNum;
		maxTime = _maxTime;
		WR = _wr;
		CR1 = _cr1;
		CR2 = _cr2;
		num_projs = _num_projs;
	}

	// ȷ���ھ�
	void define_neighbors(ProblemClass problem);

	// ��ʼ��
	void initialize(ProblemClass problem);

	// ����gBest
	void update_gBest();

	// ����ƥ�佻��
	void crossover_PMX(vector<int>& code_s1, vector<int> code_s2, int num);

	// ����������Ž�
	void get_num_sol(ProblemClass problem);


	IndividualClass Shake(ProblemClass problem, IndividualClass x, int k);

	IndividualClass Two_top_search(IndividualClass x, ProblemClass problem, int ln);

	void VNS(ProblemClass problem, int index);

	void local_search(ProblemClass problem);

	// ����������
	void Fast_DPSOClass::VNS_search(ProblemClass problem, int index);

	vector<IndividualClass> main_run(ProblemClass problem, int nicheNum);

	void resultCollation();

	bool isSameLayout(vector<vector<int>> layout1, vector<vector<int>> layout2);

};