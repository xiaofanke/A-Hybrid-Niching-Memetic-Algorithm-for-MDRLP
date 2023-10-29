#pragma once
#include "head.h"
#include "Individual.h"
#define MAXDIM   60	// max dimension
#define MAXSIZE  200    // max swarm size
#define TBSIZE 15 // Tabu list size

class Fast_DPSOClass {
public:
	int popNum; //种群大小
	int maxTime; //最大迭代时间
	double WR; //变异概率
	double CR1; //交叉概率
	double CR2; //变异概率
	int num_eval; //适应度评价次数

	vector<IndividualClass> pops; // 粒子种群
	vector<IndividualClass> allBest; // 所有最优解

	vector<IndividualClass> pBest; // 每个粒子的历史最优解
	vector<IndividualClass> gBest; // 每个小生境的粒子最优解


	vector<IndividualClass> bestSolution;
	vector<IndividualClass> allSolution;


	queue<int> tabu_list; // 禁忌表
	int tabu_index[MAXSIZE] = { 0 }; // 记录当前粒子是否被禁忌，提高查询速度
	
	int nsize; // 邻居的个数
	int num_projs; // nh 哈希函数个数
	int neighbor[MAXSIZE][MAXSIZE];

	int ls_swap[1500][2]; // 2top局部搜索所有可选的两个设备互换位置。
	vector<int> ls_index; // 所有可选的两个设备互换位置的索引

	Fast_DPSOClass(int _popNum, int _maxTime, double _wr, double _cr1, double _cr2, int _num_projs) {
		popNum = _popNum;
		maxTime = _maxTime;
		WR = _wr;
		CR1 = _cr1;
		CR2 = _cr2;
		num_projs = _num_projs;
	}

	// 确定邻居
	void define_neighbors(ProblemClass problem);

	// 初始化
	void initialize(ProblemClass problem);

	// 更新gBest
	void update_gBest();

	// 部分匹配交叉
	void crossover_PMX(vector<int>& code_s1, vector<int> code_s2, int num);

	// 获得所有最优解
	void get_num_sol(ProblemClass problem);


	IndividualClass Shake(ProblemClass problem, IndividualClass x, int k);

	IndividualClass Two_top_search(IndividualClass x, ProblemClass problem, int ln);

	void VNS(ProblemClass problem, int index);

	void local_search(ProblemClass problem);

	// 变领域搜索
	void Fast_DPSOClass::VNS_search(ProblemClass problem, int index);

	vector<IndividualClass> main_run(ProblemClass problem, int nicheNum);

	void resultCollation();

	bool isSameLayout(vector<vector<int>> layout1, vector<vector<int>> layout2);

};