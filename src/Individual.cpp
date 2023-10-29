#include "Individual.h"

bool compare_ma_dis(const vector<double>& a, const vector<double>& b) {
	return a[1] < b[1];
}


// �����ʼ��
void IndividualClass::init(ProblemClass problem) {

	vector<int> tsequence;

	// ���б���	
	for (int i = 1; i < problem.facility_num + 1; i++) {
		tsequence.push_back(i);
	}
	random_shuffle(tsequence.begin(), tsequence.end());

	// ����ϵ�
	int t = randAB(problem.facility_num / 2 - 3, problem.facility_num / 2 + 3);
	//int t = (rand() % 7) + problem.facility_num / 2 - 3;
	tsequence.insert(tsequence.begin() + t, 0);

	sequence = tsequence;	

	location.resize(problem.facility_num + 1);
	indexR.resize(problem.facility_num + 1);
}


//�����˵�λ��
void IndividualClass::adjustBreakpoint(ProblemClass problem) {

	int ln = (problem.facility_num + 1) / 2;
	int lh = node_index[0];

	// ͬ������Ҫ��������Ϊ�����������ֻ�о��ǽ�sequenceת��layout
	//// ����node_index  
	//node_index[0] = ln;
	//if (ln < hn) {
	//	for (int i = ln; i < hn; i++)
	//		node_index[sequence[i]]++;
	//}
	//else {
	//	for (int i = hn + 1; i <= ln; i++)
	//		node_index[sequence[i]]--;
	//}

	sequence.erase(sequence.begin() + lh); // ɾ��0Ԫ��
	sequence.insert(sequence.begin() + ln, 0); // ����0Ԫ��
}


// �������offset��Ŀ�꺯��ֵ
double IndividualClass::getObj_offset(int& objNum, ProblemClass problem, int num) {

	objNum++;

	// ��ȡÿ���豸��λ��
	for (int i = 0; i < problem.facility_num + 1; i++)
		node_index[sequence[i]] = i;	

	// �豸1λ�ڵڶ��У�����һ���е��豸, ȷ���豸1λ�ڵ�һ��
	if (node_index[0] < node_index[1]) { 
		vector<int> tsequence;
		for (int i = node_index[0] + 1; i < problem.facility_num + 1; i++)
			tsequence.push_back(sequence[i]);
		tsequence.push_back(0);
		for (int i = 0; i < node_index[0]; i++)
			tsequence.push_back(sequence[i]);
		sequence = tsequence;

		node_index[0] = problem.facility_num - node_index[0];
	}

	if (node_index[0] < (problem.facility_num + 1) / 2 - 3 || node_index[0] > (problem.facility_num + 1) / 2 + 3)
		adjustBreakpoint(problem);

	vector<vector<int>> tlayout;
	tlayout.resize(2);

	// �豸���ֵ���
	double tempX = 0;
	int preNode = -1;
	int tempR = 0;
	for (int i = 0; i < problem.facility_num + 1; i++) {
		int node = sequence[i];

		node_index[node] = i;

		if (node == 0) {
			preNode = -1;
			tempX = 0;
			tempR++;
			continue;
		}
		tlayout[tempR].push_back(node);

		tempX += problem.widths[node] / 2;
		if (preNode > -1) {			
			tempX += problem.minSpaces[preNode][node];
		}

		location[node] = tempX;
		indexR[node] = tempR;

		tempX += problem.widths[node] / 2;		
		preNode = node;
	}

	// ����offset	
	vector<vector<double>> interaction;

	double AllMa = 0;
	for (int i = 0; i < tlayout[0].size(); i++) {
		for (int j = 0; j < tlayout[1].size(); j++) {
			vector<double> tml = { problem.materials[tlayout[0][i]][tlayout[1][j]], location[tlayout[0][i]] - location[tlayout[1][j]] };
			interaction.push_back(tml);
			AllMa += problem.materials[tlayout[0][i]][tlayout[1][j]];
		}
	}
	sort(interaction.begin(), interaction.end(), compare_ma_dis);

	double Tma = 0;
	int indexT;
	for (int i = 0; i < interaction.size(); i++) {
		if (Tma + interaction[i][0] + 1.0e-8 >= AllMa / 2.0) {
			indexT = i;
			break;
		}
		Tma += interaction[i][0];
	}
	offset = -1 * interaction[indexT][1];


	// ����Ŀ�꺯��ֵ
	for (int i = 0; i < tlayout[0].size(); i++)
		location[tlayout[0][i]] += offset;

	//�豸1λ���豸2���Ҳ࣬��ߵ�ÿ���ϵ��豸λ��
	if (location[2] < location[1]) {
		
		// ֻ��Ҫ����sequence���У�������ֵû�б�Ҫ�����й����и��£���Ϊ��������ֻ��sequence�������
		// ������Щ����Ӱ��Ŀ�꺯���ļ��㣬���ǽ������һ��ʵ������ֵ�ø���

		reverse(sequence.begin(), sequence.begin() + node_index[0]);
		reverse(sequence.begin() + node_index[0] + 1, sequence.end());


		//// ����layout
		//reverse(tlayout[0].begin(), tlayout[0].end());
		//reverse(tlayout[1].begin(), tlayout[1].end());


		////����location
		//int node1 = tlayout[1][0]; // ����ÿ�еĵ�һ���豸��ԭ���е����һ��
		//double maxLen1 = location[node1] + problem.widths[node1] / 2.0;	

		//for (int i = 1; i < problem.facility_num + 1; i++) {
		//	location[i] = maxLen1 - location[i];
		//}

		//// ����offset
		//int node0 = tlayout[0][0];
		//double maxLen0 = location[node0] + problem.widths[node0] / 2.0;
		//offset = maxLen0 - maxLen1;
	}

	obj_offset = 0;
	for (int i = 1; i < problem.facility_num; i++) {
		for (int j = i + 1; j < problem.facility_num + 1; j++) {
			if (indexR[i] == indexR[j])
				obj_offset += fabs(location[i] - location[j]) * problem.materials[i][j];
			else
				obj_offset += (fabs(location[i] - location[j]) + problem.Cs) * problem.materials[i][j];
		}
	}

	layout = tlayout; // ���н�

	if (num > 0) {
		char saveFilename[1024];
		sprintf(saveFilename, "log/GA/obj_offset/Solution_F%d_R%d_T%d.data", problem.facility_num, problem.row_num, num);
		std::fstream fout;
		fout.open(saveFilename, std::ios::out);
		fout.precision(120);

		// Ŀ�꺯��ֵ
		fout << "optimal:" << obj_offset << endl;

		//ÿ���豸����
		fout << "sequence:" << endl;
		for (int r = 0; r < 2; r++) {
			for (int i = 0; i < layout[r].size(); i++) {
				fout << layout[r][i] - 1<< " ";
			}
			fout << endl;
		}

		// ������ֵ
		fout << "X:" << endl;
		for (int i = 1; i < problem.facility_num + 1; i++)
			fout << location[i] << " ";
		fout << endl;

		// �豸��Ӧ����
		fout << "indexR:" << endl;
		for (int i = 1; i < problem.facility_num + 1; i++) {
			fout << indexR[i] << " ";
		}
		fout << endl;
		fout.close();
	}

	return obj_offset;
}


double IndividualClass::getObj_final(ProblemClass problem, int num) {

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(env);
	try {

		// ����Լ��		
		IloNumVar min_offset = IloNumVar(env, -5000, 5000);

		// Ŀ�꺯��ֵ����
		IloNumExpr fun1(env);

		for (int i = 1; i < problem.facility_num; i++) {
			for (int j = i + 1; j < problem.facility_num + 1; j++) {
				if (indexR[i] == 0 && indexR[j] == 1)
					fun1 += (IloAbs(location[i] - location[j] + min_offset) + problem.Cs) * problem.materials[i][j];
				else if (indexR[i] == 1 && indexR[j] == 0)
					fun1 += (IloAbs(location[i] - location[j] - min_offset) + problem.Cs) * problem.materials[i][j];
				else
					fun1 += IloAbs(location[i] - location[j]) * problem.materials[i][j];
			}
		}


		model.add(IloMinimize(env, fun1));

		fun1.end();

		cplex.setParam(IloCplex::MemoryEmphasis, true);
		cplex.setParam(IloCplex::TreLim, 13000);
		cplex.setParam(IloCplex::VarSel, 3);
		cplex.setParam(IloCplex::WorkMem, 7200);
		cplex.setParam(IloCplex::NodeFileInd, 3);
		cplex.setParam(IloCplex::Threads, 1);

		cplex.setParam(IloCplex::TiLim, 40000);
		cplex.extract(model);

		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());

		IloTimer timer(env);
		timer.start();
		cplex.solve();
		timer.stop();

		offset = cplex.getValue(min_offset);
		obj_final = cplex.getObjValue();

		obj_final = round(obj_final * 1000) / 1000;

		if (num > 0) {

			if (cplex.getStatus() == IloAlgorithm::Infeasible)
				env.out() << "No Solution" << endl;
			cout << endl;

			char saveFilename[1024];
			sprintf(saveFilename, "log/GA/obj_final/Solution_F%d_R%d_T%d.data", problem.facility_num, problem.row_num, num);
			std::fstream fout;
			fout.open(saveFilename, std::ios::out);
			fout.precision(120);

			cout << "time is:" << timer.getTime() << endl;
			fout << "optimal value:" << cplex.getObjValue() << endl;

			//ÿ���豸����
			fout << "sequence:" << endl;
			for (int r = 0; r < 2; r++) {
				for (int i = 0; i < layout[r].size(); i++) {
					fout << layout[r][i] << " ";
				}
				fout << endl;
			}

			fout << "X:" << endl;
			for (int i = 1; i < problem.facility_num + 1; i++)
				if (indexR[i] == 0)
					fout << offset + location[i] << " ";
				else
					fout << location[i] << " ";
			fout << endl;

			fout << "indexR:" << endl;
			for (int i = 1; i < problem.facility_num + 1; i++) {
				fout << indexR[i] << " ";
			}
			fout << endl;
			cout << "x:" << offset << endl;
			printf("����ֵΪ%f\n", cplex.getObjValue());
			fout.close();
		}
		//return ans;
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
		cerr << ex.getMessage() << endl;
		//return ans;
	}
	catch (...) {
		cerr << "Error" << endl;
		//return ans;
	}
	cplex.end();
	model.end();
	env.end();

	return obj_final;
}