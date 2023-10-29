#include "Fast_DPSO.h"

struct objIdx {
	double obj;
	int id;
};

struct Projection
{
	double a1[MAXDIM], a2[MAXDIM];  // random Gaussian vector
	double b1, b2;          // offset
	double r1, r2;          // segment width
	double p1[MAXSIZE], p2[MAXSIZE]; // projection
	int hashid[MAXSIZE];
	map<int, vector<int> > bucket;
};
Projection proj[200];

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine gen(seed);
normal_distribution<double> dis(0, 1);

// 自定义个体比较函数
bool compare_ind(IndividualClass& a, IndividualClass& b) {
	//return a.obj_offset < b.obj_offset;
	return a.obj_final < b.obj_final;
}


// 判断两个编码串是否相同
//bool isSame(vector<int> a, vector<int> b, int num) {
//	for (int i = 0; i < num; i++) {
//		if (a[i] != b[i])
//			return 0;
//	}
//	return 1;
//}


bool compare_objIdx(objIdx& a, objIdx& b) {
	return a.obj < b.obj;
}


/**
 * 初始化种群和pBest值.
 * 
 * \param problem
 */
void Fast_DPSOClass::initialize(ProblemClass problem) {
	/* initialize parameters */
	num_eval = 0; // 目标评价次数初始化

	// 初始化种群和pBest，其中gBest在种群划分中更新
	for (int i = 0; i < popNum; i++)
	{
		IndividualClass idc = IndividualClass();
		idc.init(problem);
		idc.getObj_offset(num_eval, problem, 0);
		pops.push_back(idc);
		pBest.push_back(idc);
	}
	gBest.resize(popNum); //gBest大小初始化

	// 初始化局部搜索需要用到的数据
	int ls_num = 0;
	for (int i = 0; i < problem.facility_num; i++) {
		for (int j = i + 1; j < problem.facility_num + 1; j++) {
			ls_swap[ls_num][0] = i;
			ls_swap[ls_num][1] = j;
			ls_index.push_back(ls_num);
			ls_num++;
		}
	}
}


// 生成邻居
void Fast_DPSOClass::define_neighbors(ProblemClass problem) {

	// 构造哈希函数
	for (int i = 0; i < num_projs; i++) {

		// 构造两个高斯分布向量
		for (int j = 0; j < problem.facility_num + 1; j++) {
			proj[i].a1[j] = dis(gen);
			proj[i].a2[j] = dis(gen);
		}

		// 初始化索引(注意，下面这三行好像没有意义，属于多余的行)
		map<int, vector<int> >::iterator iter;
		for (iter = proj[i].bucket.begin(); iter != proj[i].bucket.end(); iter++)
			iter->second.clear();

		// 初始化索引
		proj[i].bucket.clear();

		// XV的上下界
		double min1, min2, max1, max2;
		min1 = min2 = 1e300;
		max1 = max2 = -1e300;

		// 计算V*X， 其中X是进行了归一化，即将当前值减去了(Xmax-Xmin)/2
		for (int j = 0; j < popNum; j++) {
			proj[i].p1[j] = 0;
			proj[i].p2[j] = 0;
			for (int k = 0; k < problem.facility_num + 1; k++) {
				proj[i].p1[j] += proj[i].a1[k] * (pBest[j].sequence[k] - (problem.facility_num + 1) / 2.0); // 注意这里使用的是Pbest，不是当前例子
				proj[i].p2[j] += proj[i].a2[k] * (pBest[j].sequence[k] - (problem.facility_num + 1) / 2.0);
			}
			if (proj[i].p1[j] > max1) max1 = proj[i].p1[j];
			if (proj[i].p1[j] < min1) min1 = proj[i].p1[j];

			if (proj[i].p2[j] > max2) max2 = proj[i].p2[j];
			if (proj[i].p2[j] < min2) min2 = proj[i].p2[j];
		}

		double d = 20;  // number of buckets， nb
		proj[i].r1 = (max1 - min1) / (d);
		proj[i].r2 = (max2 - min2) / (d);
		proj[i].b1 = myRand() * proj[i].r1;
		proj[i].b2 = myRand() * proj[i].r2;

		// 公式7，Z的值定义为50
		for (int j = 0; j < popNum; j++) {
			double id = ceil((proj[i].p1[j] + proj[i].b1) / proj[i].r1) * 50 + ceil((proj[i].p2[j] + proj[i].b2) / proj[i].r2);
			proj[i].hashid[j] = id;
			proj[i].bucket[id].push_back(j);
		}
	}


	// 小生境生成
	for (int i = 0; i < popNum; i++) {
		for (int j = 0; j < nsize; j++) {
			int np = rand() % num_projs;
			int id = proj[np].hashid[i];
			int bs = proj[np].bucket[id].size();
			neighbor[i][j] = proj[np].bucket[id][rand() % bs];
		}
	}
}


// 更新gBest
void Fast_DPSOClass::update_gBest() {
	for (int i = 0; i < popNum; i++) {
		//gBest[i] = pBest[i]; // 设备i的邻居应该也包含自己
		for (int j = 0; j < nsize; j++) {
			int neId = neighbor[i][j];
			if (j == 0) {
				gBest[i] = pBest[neId];
				continue;
			}			
			if (pBest[neId].obj_offset < gBest[i].obj_offset) {
				gBest[i] = pBest[neId];
			}
		}
	}
}


//PMX交叉操作 部分匹配交叉
void Fast_DPSOClass::crossover_PMX(vector<int>& code_s1, vector<int> code_s2, int num) {
	
	//随机生成两个交叉点，保证lw小于hg
	int lw = rand() % num;
	int hg = rand() % num;

	if (lw > hg) {
		lw = lw ^ hg;
		hg = lw ^ hg;
		lw = lw ^ hg;
	}
	hg++; // 包含hg[lw, hg]

	map<int, int> map2T1;

	//替换操作	
	for (int i = lw; i < hg; i++) {		
		map2T1.insert(pair<int, int>(code_s2[i], code_s1[i]));
		code_s1[i] = code_s2[i];		
	}

	map<int, int>::iterator it;
	//整理替换解
	for (int i = 0; i < lw; i++) {
		int tn = code_s1[i];
		it = map2T1.find(tn);
		while (it != map2T1.end()) {
			tn = it->second;
			it = map2T1.find(tn);
		}
		code_s1[i] = tn;		
	}

	for (int i = hg; i < num; i++) {
		int tn = code_s1[i];
		it = map2T1.find(tn);
		while (it != map2T1.end()) {
			tn = it->second;
			it = map2T1.find(tn);
		}
		code_s1[i] = tn;		
	}
}


void Fast_DPSOClass::get_num_sol(ProblemClass problem) {
	sort(pBest.begin(), pBest.end(), compare_ind);
	allBest.push_back(pBest[0]);

	for (int i = 1; i < popNum; i++) {
		if (fabs(pBest[0].obj_offset - pBest[i].obj_offset) < 1.0e-2) {
			pBest[i].getObj_offset(num_eval, problem, 0);

			int flg = 1;
			for (int j = 0; j < allBest.size(); j++) {
				if (isSameLayout(pBest[i].layout, allBest[j].layout) == 1) {
					flg = 0;
					break;
				}
			}
			if (flg == 1) {				
				allBest.push_back(pBest[i]);
			}
		}
		else
			break;
	}
}


IndividualClass Fast_DPSOClass::Two_top_search(IndividualClass x, ProblemClass problem, int ln) {

	int k = 1;
	while (k) {
		k = 0;
		random_shuffle(ls_index.begin(), ls_index.end());

		//cout << "-----1-----" << endl;

		for (int i = 0; i < problem.facility_num * 2; i++) {
			int rd = rand() % 2;
			IndividualClass x1 = x;
			if (rd == 1) {
				int a = ls_swap[ls_index[i]][0];
				int b = ls_swap[ls_index[i]][1];

				//cout << "-----2-----" << endl;

				x1.sequence[a] = x1.sequence[a] ^ x1.sequence[b];
				x1.sequence[b] = x1.sequence[a] ^ x1.sequence[b];
				x1.sequence[a] = x1.sequence[a] ^ x1.sequence[b];

				//cout << "-----2.1-----" << endl;
				/*for (int j = 0; j < x1.sequence.size(); j++)
					cout << x1.sequence[j] << " ";
				cout << endl;*/

				x1.getObj_offset(num_eval, problem, 0);

				//cout << "-----2.2-----" << endl;

				if (x1.obj_offset < x.obj_offset) {
					x = x1;
					k = 1;
					break;
				}
			}
			else {

				//cout << "-----3-----" << endl;

				int a = ls_swap[ls_index[i]][0];
				int b = ls_swap[ls_index[i]][1];


				int temp = x1.sequence[a];
				x1.sequence.erase(x1.sequence.begin() + a);
				x1.sequence.insert(x1.sequence.begin() + b, temp);

				x1.getObj_offset(num_eval, problem, 0);

				if (x1.obj_offset < x.obj_offset) {
					x = x1;
					k = 1;
					break;
				}
			}
		}
	}
	return x;
}


IndividualClass Fast_DPSOClass::Shake(ProblemClass problem, IndividualClass x, int k) {
	if (k == 0) { // 互换两个设备
		int a = rand() % (problem.facility_num + 1);
		int b = rand() % (problem.facility_num + 1);

		while (a == b) {
			b = rand() % (problem.facility_num + 1);
		}

		x.sequence[a] = x.sequence[a] ^ x.sequence[b];
		x.sequence[b] = x.sequence[a] ^ x.sequence[b];
		x.sequence[a] = x.sequence[a] ^ x.sequence[b];
		
		x.getObj_offset(num_eval, problem, 0);

		return x;
	}
	else if (k == 1) { // 选择一个设备放到一个可行位置
		int a = rand() % (problem.facility_num + 1);
		int b = rand() % (problem.facility_num + 1);

		while (a == b) {
			b = rand() % (problem.facility_num + 1);
		}

		int temp = x.sequence[a];
		x.sequence.erase(x.sequence.begin() + a);
		x.sequence.insert(x.sequence.begin() + b, temp);

		x.getObj_offset(num_eval, problem, 0);

		return x;
	}
	else if (k == 2) { // 互换三个设备
		int a = rand() % (problem.facility_num + 1);
		int b = rand() % (problem.facility_num + 1);
		int c = rand() % (problem.facility_num + 1);

		while (a == b) {
			b = rand() % (problem.facility_num + 1);
		}
		while (a == c || b == c) {
			c = rand() % (problem.facility_num + 1);
		}

		int temp = x.sequence[a];
		x.sequence[a] = x.sequence[b];
		x.sequence[b] = x.sequence[c];
		x.sequence[c] = temp;

		x.getObj_offset(num_eval, problem, 0);

		return x;
	}
	else if (k == 3) { // 随机选择一个子序列颠倒位置

		int n = x.sequence.size();

		int i = rand() % n;
		int q = rand() % (n / 4 - n / 8) + (n / 8 + 1);
		if (n < 20)
			q = rand() % 2 + 3;
		int j = i + q - 1;
		if (j >= n)
			j = j - n;

		for (int k = 1; k < q / 2; k++) {
			int temp = x.sequence[i];
			x.sequence[i] = x.sequence[j];
			x.sequence[j] = temp;

			i++;
			if (i >= n)
				i = 0;
			j--;
			if (j < 0)
				j = n - 1;
		}

		x.getObj_offset(num_eval, problem, 0);

		return x;
	}
}


void Fast_DPSOClass::VNS(ProblemClass problem, int index) {

	double rd = myRand();

	IndividualClass x = pBest[index];
	if (rd < 0.5)
		x = pops[index];

	int k = 0;
	while (k <= 3) {
		IndividualClass x1 = Shake(problem, x, k);
		IndividualClass x2 = Two_top_search(x1, problem, problem.facility_num);

		if (x2.obj_offset < x.obj_offset) {
			x = x2;
			k = 1;
		}
		else {
			k++;
		}
	}

	if (x.obj_offset < pBest[index].obj_offset) {
		pBest[index] = x;
	}
	if (rd < 0.5 && x.obj_offset < pops[index].obj_offset) {
		pops[index] = x;
	}

}


void Fast_DPSOClass::local_search(ProblemClass problem) {
	vector<objIdx> tempObj;
	for (int i = 0; i < popNum; i++) {
		objIdx idx = { pBest[i].obj_offset, i };
		tempObj.push_back(idx);
	}
	sort(tempObj.begin(), tempObj.end(), compare_objIdx);

	//先调整pbest，去掉重复解
	double tempF = tempObj[0].obj; //当前最小值
	//int tempId = tempObj[0].id;
	int index = 0; // 当前最小值的开始索引
	for (int i = 1; i < popNum; i++) {

		if (fabs(tempF - pBest[tempObj[i].id].obj_offset) < 1.0e-5) { //布局成本等于当前最小值
			int flg = 0;
			for (int j = index; j < i; j++) {
				if (isSameLayout(pBest[j].layout, pBest[tempObj[i].id].layout) == 1) { // 判断布局结构是否相同
					flg = 1;
					break;
				}
			}
			if (flg) {
				IndividualClass idc = IndividualClass();
				idc.init(problem);
				idc.getObj_offset(num_eval, problem, 0);

				pops[tempObj[i].id] = idc;
				pBest[tempObj[i].id] = idc;
			}
		}
		else { // 更换最小值和最小值索引
			tempF = tempObj[i].obj;
			//tempId = tempObj[i].id;
			index = i;
		}
	}

	tempObj.clear();
	for (int i = 0; i < popNum; i++) {
		objIdx idx = { pBest[i].obj_offset, i };
		tempObj.push_back(idx);
	}
	sort(tempObj.begin(), tempObj.end(), compare_objIdx);


	int Num_local = 0;
	for (int i = 0; i < popNum; i++) {
		if (Num_local >= 5)
			break;
		if (tabu_index[tempObj[i].id] == 1)
			continue;

		// 修改禁忌表
		tabu_list.push(tempObj[i].id);
		tabu_index[tempObj[i].id] = 1;

		if (tabu_list.size() > TBSIZE) // 禁忌表的长度超过阈值，删除栈顶元素，修改tabu_index
		{
			int tempIdx = tabu_list.front();
			tabu_index[tempIdx] = 0;
			tabu_list.pop();
		}

		VNS(problem, tempObj[i].id);
		//swap_search(problem, tempObj[i].id);
		Num_local++;
	}
}


vector<IndividualClass> Fast_DPSOClass::main_run(ProblemClass problem, int ns) {

	// 初始化
	initialize(problem);
	
	// 迭代优化
	for (int iter = 0; iter < maxTime; iter++) {

		//cout << iter << " " << num_eval << endl;

		// 更新nsize
		if (iter <= 0.25 * maxTime)
			nsize = 2;
		else if (iter <= 0.5 * maxTime)
			nsize = 2 + (ns-2) / 3;
		else if (iter <= 0.75 * maxTime)
			nsize = 2 + (ns - 2) / 3 * 2;
		else
			nsize = ns;

		/*if (iter % (2 * num_projs) == 0)
			define_neighbors(problem);*/

		if (iter % (num_projs) == 0)
			define_neighbors(problem);

		update_gBest(); // 更新粒子的gBest

		// 粒子更新
		for (int i = 0; i < popNum; i++) {
			//当前粒子进行突变
			if (myRand() < WR) {
				// 随机选择一个设备位置，2-facility_num
				int a = rand() % (problem.facility_num + 1);
				int b = rand() % (problem.facility_num + 1);
				while (a == b) {
					b = rand() % (problem.facility_num + 1);
				}

				// 互换a，b位置上的设备
				pops[i].sequence[a] = pops[i].sequence[a] ^ pops[i].sequence[b];
				pops[i].sequence[b] = pops[i].sequence[a] ^ pops[i].sequence[b];
				pops[i].sequence[a] = pops[i].sequence[a] ^ pops[i].sequence[b];
			}

			//变异粒子与pBest交叉
			if (myRand() < CR1) {
				crossover_PMX(pops[i].sequence, pBest[i].sequence, problem.facility_num + 1);
			}

			//交叉后粒子与gBest交叉
			if (myRand() < CR2) {
				crossover_PMX(pops[i].sequence, gBest[i].sequence, problem.facility_num + 1);
			}

			// 计算目标函数值
			pops[i].getObj_offset(num_eval, problem, 0);

			// 更新粒子的历史最优解
			if (pops[i].obj_offset < pBest[i].obj_offset) {
				pBest[i] = pops[i];
				tabu_index[i] = 0; // 移除禁忌表
			}
		}
		local_search(problem);
	}

	for (int i = 0; i < pBest.size(); i++)
		pBest[i].getObj_final(problem, 0);

	//// 统计所有最优解
	//cout << "目标评价次数：" << num_eval << endl;
	get_num_sol(problem);

	return allBest;
}


/**
 * 结果统计，包括最优解（含统计），所有解（含统计）.
 * 
 */
void Fast_DPSOClass::resultCollation() {

	vector<IndividualClass> bs; // 最优解
	vector<IndividualClass> as; // 所有解


	sort(pBest.begin(), pBest.end(), compare_ind);

	double minF = pBest[0].obj_final;
	bs.push_back(pBest[0]);
	as.push_back(pBest[0]);

	for (int i = 1; i < pBest.size(); i++) {
		if (fabs(pBest[i].obj_final - minF) < 1.0e-5) { // 与最优解的目标函数值相同
			int flg = 1;
			for (int j = 0; j < bs.size(); j++) {
				if (isSameLayout(bs[j].layout, pBest[i].layout)) {
					flg = 0;
					break;
				}				
			}
			if (flg) {
				bs.push_back(pBest[i]);
				as.push_back(pBest[i]);
			}
		}
		else {
			int flg = 1;
			for (int j = 0; j < as.size(); j++) {
				if (isSameLayout(as[j].layout, pBest[i].layout)) {
					flg = 0;
					break;
				}
			}
			if (flg) {
				as.push_back(pBest[i]);
			}
		}
	}

	bestSolution = bs;
	allSolution = as;

}


/**
 * 判断两个解（布局）是否相同，主要考虑两行的设备排列.
 *
 * \param layout1 布局1
 * \param layout2 布局2
 * \return
 */
bool Fast_DPSOClass::isSameLayout(vector<vector<int>> layout1, vector<vector<int>> layout2) {

	//判断完全一致
	//  1 2 3 4    1 2 3 4
	//  5 6 7 8    5 6 7 8
	if (layout1[0].size() == layout2[0].size()) {

		int tf = 1;
		for (int i = 0; i < layout1[0].size(); i++) // 第一行对比
			if (layout1[0][i] != layout2[0][i]) {
				tf = 0;
				break;
			}
		if (tf) // 第二行对比
			for (int i = 0; i < layout1[1].size(); i++)
				if (layout1[1][i] != layout2[1][i]) {
					tf = 0;
					break;
				}
		if (tf)
			return true;
	}

	//判断对称互换
	//  1 2 3 4    4 3 2 1
	//  5 6 7 8    8 7 6 5
	if (layout1[0].size() == layout2[0].size()) {

		int tf = 1;
		for (int i = 0; i < layout1[0].size(); i++) // 第一行对比
			if (layout1[0][i] != layout2[0][layout2[0].size() - 1 - i]) {
				tf = 0;
				break;
			}
		if (tf) // 第二行对比
			for (int i = 0; i < layout1[1].size(); i++)
				if (layout1[1][i] != layout2[1][layout2[1].size() - 1 - i]) {
					tf = 0;
					break;
				}
		if (tf)
			return true;
	}

	//判断两行互换
	//  1 2 3 4    5 6 7 8
	//  5 6 7 8    1 2 3 4
	if (layout1[0].size() == layout2[1].size()) {

		int tf = 1;
		for (int i = 0; i < layout1[0].size(); i++) // 第一行对比
			if (layout1[0][i] != layout2[1][i]) {
				tf = 0;
				break;
			}
		if (tf) // 第二行对比
			for (int i = 0; i < layout1[1].size(); i++)
				if (layout1[1][i] != layout2[0][i]) {
					tf = 0;
					break;
				}
		if (tf)
			return true;
	}

	//判断对称且两行互换
	//  1 2 3 4    8 7 6 5
	//  5 6 7 8    4 3 2 1
	if (layout1[0].size() == layout2[1].size()) {

		int tf = 1;
		for (int i = 0; i < layout1[0].size(); i++) // 第一行对比
			if (layout1[0][i] != layout2[1][layout2[1].size() - 1 - i]) {
				tf = 0;
				break;
			}
		if (tf) // 第二行对比
			for (int i = 0; i < layout1[1].size(); i++)
				if (layout1[1][i] != layout2[0][layout2[0].size() - 1 - i]) {
					tf = 0;
					break;
				}
		if (tf)
			return true;
	}

	return false;
}
