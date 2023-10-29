#include "head.h"
#include "Fast_DPSO.h"

int main() {
	srand((unsigned)time(NULL));
	clock_t startTime, endTime;

	//vector<vector<int>> paprameters = { // 参数1
		/*{8, 2, 100},
		{ 8, 3, 100 },
		{ 10, 2, 300 },
		{ 10, 3, 300 },
		{ 12, 3, 700 },
		{ 12, 4, 700 },
		{ 16, 3, 1000 },
		{ 16, 4, 1000 },
		{ 18, 4, 1500 },
		{ 18, 5, 1500 },
		{ 20, 5, 2000 },
		{ 20, 6, 2000 },
		{ 26, 6, 2500 },
		{ 30, 6, 3000 },*/
	//};

	vector<vector<int>> paprameters = { // 参数2
		{8, 2, 100},
		{8, 3, 100},
		{10, 2, 300},
		{10, 3, 300},
		{12, 3, 700},
		{12, 4, 700},
		{16, 3, 1000},
		{16, 4, 1000},
		{18, 4, 1500},
		{18, 5, 1500},	
		{20, 5, 2000},
		{20, 6, 2000},
		{26, 6, 2500},
		{30, 6, 3000},			
	};

	

	char fileDir[1024];
	char fileDir_bestAns_solutions[1024];
	char fileDir_bestAns_statistic[1024];

	char fileDir_allAns_solutions[1024];
	char fileDir_allAns_statistic[1024];

	int nn = 24;

	for (int nbs = 50; nbs <= 50; nbs += 5) {

		//for (int p = 9; p <= 9; p++){
		for (int p = 0; p < paprameters.size(); p++) {
			//for (int p = nn; p <= nn; p++) {

			sprintf(fileDir, "data_multi_sol/facility_F%d_R2_S%d.txt", paprameters[p][0], paprameters[p][1]);
			sprintf(fileDir_bestAns_statistic, "solution/best/facility_F%d_R2_S%d_bestStatistic.txt", paprameters[p][0], paprameters[p][1]);
			sprintf(fileDir_bestAns_solutions, "solution/best/facility_F%d_R2_S%d_beestSolutions.txt", paprameters[p][0], paprameters[p][1]);

			sprintf(fileDir_allAns_statistic, "solution/all/facility_F%d_R2_S%d_allStatistic.txt", paprameters[p][0], paprameters[p][1]);
			sprintf(fileDir_allAns_solutions, "solution/all/facility_F%d_R2_S%d_allSolutions.txt", paprameters[p][0], paprameters[p][1]);

			/*sprintf(fileDir_bestAns_statistic, "solution/best/facility_F%d_R2_S%d_bestStatistic_nbs_%d.txt", paprameters[p][0], paprameters[p][1], nbs);
			sprintf(fileDir_bestAns_solutions, "solution/best/facility_F%d_R2_S%d_beestSolutions_nbs_%d.txt", paprameters[p][0], paprameters[p][1], nbs);

			sprintf(fileDir_allAns_statistic, "solution/all/facility_F%d_R2_S%d_allStatistic_nbs_%d.txt", paprameters[p][0], paprameters[p][1], nbs);
			sprintf(fileDir_allAns_solutions, "solution/all/facility_F%d_R2_S%d_allSolutions_nbs_%d.txt", paprameters[p][0], paprameters[p][1], nbs);*/


			ProblemClass problem = ProblemClass(fileDir);
						
			
			for (int IN = 0; IN < 20; IN++) {

				std::fstream bestFout; // 最优解
				std::fstream bestStaticFout; // 最优解统计
				std::fstream allFout; // 所有解
				std::fstream allStaticFout; // 所有解统计

				// 第一次覆盖运行覆盖写入
				if (IN == 0) {
					bestFout.open(fileDir_bestAns_solutions, std::ios::out);
					bestStaticFout.open(fileDir_bestAns_statistic, std::ios::out);
					allFout.open(fileDir_allAns_solutions, std::ios::out);
					allStaticFout.open(fileDir_allAns_statistic, std::ios::out);
				}
				else { // 随后追加写入
					bestFout.open(fileDir_bestAns_solutions, std::ios::app);
					bestStaticFout.open(fileDir_bestAns_statistic, std::ios::app);
					allFout.open(fileDir_allAns_solutions, std::ios::app);
					allStaticFout.open(fileDir_allAns_statistic, std::ios::app);
				}

				// 书写格式定义
				bestFout << fixed << setprecision(5);
				bestStaticFout << fixed << setprecision(5);
				allFout << fixed << setprecision(5);
				allStaticFout << fixed << setprecision(5);

				if (IN == 0) {
					bestFout << "--------facility_" << paprameters[p][0] << "soution_" << paprameters[p][1] << "--------" << endl;
					bestStaticFout << "--------facility_" << paprameters[p][0] << "soution_" << paprameters[p][1] << "--------" << endl;
					allFout << "--------facility_" << paprameters[p][0] << "soution_" << paprameters[p][1] << "--------" << endl;
					allStaticFout << "--------facility_" << paprameters[p][0] << "soution_" << paprameters[p][1] << "--------" << endl;
				}

			
				bestFout << "第" << IN + 1 << "次运行：" << endl;
				allFout << "第" << IN + 1 << "次运行：" << endl;

				Fast_DPSOClass fdp = Fast_DPSOClass(100, paprameters[p][2], 0.2, 0.6, 0.8, 10);

				startTime = clock();
				fdp.main_run(problem, nbs);
				endTime = clock();

				fdp.resultCollation();

				cout << "时间" << double(endTime - startTime) / CLOCKS_PER_SEC << endl;

				bestStaticFout << "独立运行次数: " << IN + 1 << " ";
				bestStaticFout << "解的个数：" << fdp.bestSolution.size() << " ";
				bestStaticFout << "评价次数：" << fdp.num_eval << " ";
				bestStaticFout << "运行时间：" << double(endTime - startTime) / CLOCKS_PER_SEC << " ";
				bestStaticFout << "所有解的目标值：" << fdp.bestSolution[0].obj_final << endl;

				allStaticFout << "独立运行次数: " << IN + 1 << " ";
				int AN = pow(2, paprameters[p][1] - 1);
				if (fdp.allSolution.size() < AN)
					AN = fdp.allSolution.size();
				allStaticFout << "解的个数：" << AN << " ";
				allStaticFout << "评价次数：" << fdp.num_eval << " ";
				allStaticFout << "运行时间：" << double(endTime - startTime) / CLOCKS_PER_SEC << " ";
				allStaticFout << "所有解的目标值：";

				bestFout << "解的个数：" << fdp.bestSolution.size() << " " << endl;
				bestFout << "评价次数：" << fdp.num_eval << " " << endl;;
				bestFout << "运行时间：" << double(endTime - startTime) / CLOCKS_PER_SEC << " " << endl;
				bestFout << "目标函数值：" << fdp.bestSolution[0].obj_final << endl;

				allFout << "解的个数：" << fdp.bestSolution.size() << " " << endl;
				allFout << "评价次数：" << fdp.num_eval << " " << endl;;
				allFout << "运行时间：" << double(endTime - startTime) / CLOCKS_PER_SEC << " " << endl;
				allFout << "目标函数值：" << fdp.bestSolution[0].obj_final << endl;

				for (int f = 0; f < fdp.bestSolution.size(); f++) {
					cout << "目标值：" << fdp.bestSolution[0].obj_final << "  评价次数：" << fdp.num_eval << endl;

					bestFout << "sequence:" << endl;
					for (int i = 0; i < 2; i++) {
						for (int j = 0; j < fdp.bestSolution[f].layout[i].size(); j++)
							bestFout << fdp.bestSolution[f].layout[i][j] << " ";
						bestFout << endl;
					}

					bestFout << "X:" << endl;
					for (int i = 1; i < problem.facility_num + 1; i++)
						bestFout << fdp.bestSolution[f].location[i] << " ";
					bestFout << endl;

					bestFout << "indexR:" << endl;
					for (int i = 1; i < problem.facility_num + 1; i++)
						bestFout << fdp.bestSolution[f].indexR[i] << " ";
					bestFout << endl;
					bestFout << endl;
				}

				for (int f = 0; f < AN; f++) {
					//cout << "目标值：" << allPs[0].fitness << "  评价次数：" << num_eval << endl;

					allStaticFout << fdp.allSolution[f].obj_final << " ";

					allFout << "sequence:" << endl;
					for (int i = 0; i < 2; i++) {
						for (int j = 0; j < fdp.allSolution[f].layout[i].size(); j++)
							allFout << fdp.allSolution[f].layout[i][j] << " ";
						allFout << endl;
					}

					allFout << "X:" << endl;
					for (int i = 1; i < problem.facility_num + 1; i++)
						allFout << fdp.allSolution[f].location[i] << " ";
					allFout << endl;

					allFout << "indexR:" << endl;
					for (int i = 1; i < problem.facility_num + 1; i++)
						allFout << fdp.allSolution[f].indexR[i] << " ";
					allFout << endl;
					allFout << endl;
				}

				allFout << endl;
				allStaticFout << endl;
				cout << "---------------" << endl;
				cout << endl;

				bestFout.close();
				bestStaticFout.close();
				allFout.close();
				allStaticFout.close();
			}			
		}

	}
	
	system("pause");
	return 0;
}

