/*********************************************************************************

							 Triangluar Face Check

							 Updating in 2021/05/21

							   By Dr. Chenlei Lv

			The functions includes:
			1. Check the triangular face connection;
			   1.1 vertex in an edge
			   1.2 (one edge, shared by three face)a face is inside an other face

**********************************************************************************/
#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm> 
using namespace std;
#pragma once

class Point_In_Edge {
public:
	int b1;
	int b2;
	int index;
};


class Triangluar_Check {

private:

	vector<vector<double>> pointSet;
	vector<vector<int>> faceSet;
	vector<vector<int>> pointNeighbor;//store the triangular face 
	vector<vector<int>> pointEdge;//store the edge	
	vector<int> errorEdge;//record the edge with repeat subedge
	vector<bool> errorEdgeJ;//record the error edge (processed or not) 

public:

	void TC_init(vector<vector<double>> pointSet_input, vector<vector<int>> faceSet_input) {

		cout << "Triangluar_Check init run!" << endl;
		pointSet = pointSet_input;
		faceSet = faceSet_input;

		cout << "Point:" << pointSet.size() << endl;
		cout << "Face:" << faceSet.size() << endl;

		pointNeighbor.resize(pointSet.size());
		for (int i = 0; i < faceSet.size(); i++) {
			int b1 = faceSet[i][0];
			int b2 = faceSet[i][1];
			int b3 = faceSet[i][2];
			//cout << b1 <<","<< b2 << "," << b3 << endl;
			pointNeighbor[b1].push_back(b2);
			pointNeighbor[b1].push_back(b3);
			pointNeighbor[b2].push_back(b3);
			pointNeighbor[b2].push_back(b1);
			pointNeighbor[b3].push_back(b1);
			pointNeighbor[b3].push_back(b2);
		}

		TC_Removerepeat_Neighbor();
		TC_init_Edge();

	}

	void TC_start(string objname, bool saveJudge) {

		cout << "Triangular_Check start." << endl;

		//TC_Print();

		TC_check_Point_In_Edge();

		//second round check, avoid an edge share more than two points
		bool judge = true;
		while (judge) {
			judge = TC_check_Point_In_Edge_doubleCheck();
		}

		//TC_Print();

		TC_Neighbor2Face();

		//if (!saveJudge) {
			//string obj_TC = "Data\\subRemesh\\" + objname + "_triangularCheck.obj";
			//TC_SaveOBJ(obj_TC, pointSet, faceSet);
		//}
		cout << endl;
		cout << "Triangular_Check finish." << endl;

	}

	vector<vector<int>> TC_Get_faceSet() {

		return faceSet;

	}

	vector<vector<double>> TC_Get_pointSet() {

		return pointSet;

	}

private:

	//*************************************init*****************************************
	void TC_Removerepeat_Neighbor() {

		for (int i = 0; i < pointNeighbor.size(); i++) {

			vector<int> pointNeighbor_i;
			for (int j = 0; j < pointNeighbor[i].size() / 2; j++) {
				int b2 = pointNeighbor[i][2 * j];
				int b3 = pointNeighbor[i][2 * j + 1];
				if (TC_exist_Neighbor(b2, b3, pointNeighbor_i)) {
					continue;
				}
				else {
					pointNeighbor_i.push_back(b2);
					pointNeighbor_i.push_back(b3);
				}
			}
			pointNeighbor[i].clear();
			pointNeighbor[i] = pointNeighbor_i;

		}

	}

	bool TC_exist_Neighbor(int b2, int b3, vector<int> pointNeighbor_i) {

		for (int i = 0; i < pointNeighbor_i.size() / 2; i++) {

			int b2i = pointNeighbor_i[2 * i];
			int b3i = pointNeighbor_i[2 * i + 1];
			if (b2 == b2i && b3 == b3i) {
				return true;
			}
			else if (b2 == b3i && b3 == b2i) {
				return true;
			}
			else {
				continue;
			}

		}
		return false;

	}

	void TC_check_Point_In_Edge() {

		vector<Point_In_Edge> checkEdge;
		for (int i = 0; i < pointEdge.size(); i++) {

			if (pointEdge.size() - i % 1000 == 0) {
				cout << pointEdge.size() - i << ";";
			}

			int b1 = i;
			vector<int> pointEdge_i = pointEdge[i];

			for (int j = 0; j < pointEdge_i.size(); j++) {
				int b2 = pointEdge_i[j];
				if (b2 < b1) {
					continue;
				}
				else {
					vector<Point_In_Edge> checkEdge_ij = TV_check_Edge(b1, b2);

					if (checkEdge_ij.size() > 0) {
						checkEdge.insert(checkEdge.end(), checkEdge_ij.begin(), checkEdge_ij.end());
					}
				}
			}
		}
		//record the edge with repeat subedge
		TC_Update_erroredge(checkEdge);
		cout << "Init:" << checkEdge.size();
		//errorEdgeJ.resize(errorEdge.size(), false);
		TV_reconnect_Error_Edge(checkEdge);
	}

	void TC_Update_erroredge(vector<Point_In_Edge> checkEdge) {

		vector<int> checkPointList;
		for (int i = 0; i < checkEdge.size(); i++) {
			int b1_i = checkEdge[i].b1;
			int b2_i = checkEdge[i].b2;
			int index_i = checkEdge[i].index;

			vector<int> b1N = pointEdge[b1_i];
			vector<int> b2N = pointEdge[b2_i];
			vector<int> indexN = pointEdge[index_i];

			if (!TC_exist_Edge(b1_i, errorEdge)) {
				errorEdge.push_back(b1_i);
			}
			if (!TC_exist_Edge(b2_i, errorEdge)) {
				errorEdge.push_back(b2_i);
			}
			if (!TC_exist_Edge(index_i, errorEdge)) {
				errorEdge.push_back(index_i);
			}

			for (int j = 0; j < b1N.size(); j++) {
				if (!TC_exist_Edge(b1N[j], errorEdge)) {
					errorEdge.push_back(b1N[j]);
				}

			}
			for (int j = 0; j < b2N.size(); j++) {
				if (!TC_exist_Edge(b2N[j], errorEdge)) {
					errorEdge.push_back(b2N[j]);
				}

			}
			for (int j = 0; j < indexN.size(); j++) {
				if (!TC_exist_Edge(indexN[j], errorEdge)) {
					errorEdge.push_back(indexN[j]);
				}
			}
		}
	}

	bool TC_check_Point_In_Edge_doubleCheck() {

		pointEdge.clear();
		TC_init_Edge();

		vector<Point_In_Edge> checkEdge;
		vector<int> checkPointList = errorEdge;

		sort(checkPointList.begin(), checkPointList.end());

		for (int i = 0; i < checkPointList.size(); i++) {
			int b1 = checkPointList[i];
			vector<int> pointEdge_i = pointEdge[b1];
			for (int j = 0; j < pointEdge_i.size(); j++) {
				int b2 = pointEdge_i[j];
				if (b2 < b1) {
					continue;
				}
				else {
					vector<Point_In_Edge> checkEdge_ij = TV_check_Edge(b1, b2);
					if (checkEdge_ij.size() > 0) {
						checkEdge.insert(checkEdge.end(), checkEdge_ij.begin(), checkEdge_ij.end());
					}
				}
			}
		}

		//errorEdgeJ.resize(errorEdge.size(), false);
		if (checkEdge.size() == 0) {
			return false;
		}
		else {
			cout << "doubleCheck:" << checkEdge.size() << ";";
			TC_Update_erroredge(checkEdge);
			TV_reconnect_Error_Edge(checkEdge);
			return true;
		}
	}

	void TV_reconnect_Error_Edge(vector<Point_In_Edge> checkEdge) {

		for (int i = 0; i < checkEdge.size(); i++) {

			Point_In_Edge checkEdge_i = checkEdge[i];
			bool judge_i = TV_reconnect_Error_Edge_singleFace(checkEdge_i);
			//errorEdgeJ[i] = judge_i;

		}

	}

	bool TV_reconnect_Error_Edge_singleFace(Point_In_Edge pie) {

		int b1 = pie.b1;
		int b2 = pie.b2;
		int b3 = -1;
		int index = pie.index;

		for (int i = 0; i < pointNeighbor[b1].size() / 2; i++) {
			int b1n1 = pointNeighbor[b1][2 * i];
			int b1n2 = pointNeighbor[b1][2 * i + 1];
			if (b1n1 == b2) {
				b3 = b1n2;
			}
			else if (b1n2 == b2) {
				b3 = b1n1;
			}
			else {
				continue;
			}
		}

		if (b3 == index) {
			return false;
		}

		if (b3 != -1) {
			//update pointNeighbor b1
			vector<int> pointNeighbor_b1;
			for (int i = 0; i < pointNeighbor[b1].size() / 2; i++) {
				int b11 = pointNeighbor[b1][2 * i];
				int b12 = pointNeighbor[b1][2 * i + 1];
				if (b11 == b2 && b12 == b3) {
					pointNeighbor_b1.push_back(index);
					pointNeighbor_b1.push_back(b3);
				}
				else if (b12 == b2 && b11 == b3) {
					pointNeighbor_b1.push_back(b3);
					pointNeighbor_b1.push_back(index);
				}
				else {
					pointNeighbor_b1.push_back(b11);
					pointNeighbor_b1.push_back(b12);
				}
			}
			pointNeighbor[b1].clear();
			pointNeighbor[b1] = pointNeighbor_b1;

			//update pointNeighbor b2
			vector<int> pointNeighbor_b2;
			for (int i = 0; i < pointNeighbor[b2].size() / 2; i++) {
				int b21 = pointNeighbor[b2][2 * i];
				int b22 = pointNeighbor[b2][2 * i + 1];
				if (b21 == b1 && b22 == b3) {
					pointNeighbor_b2.push_back(index);
					pointNeighbor_b2.push_back(b3);
				}
				else if (b21 == b3 && b22 == b1) {
					pointNeighbor_b2.push_back(b3);
					pointNeighbor_b2.push_back(index);
				}
				else {
					pointNeighbor_b2.push_back(b21);
					pointNeighbor_b2.push_back(b22);
				}
			}
			pointNeighbor[b2].clear();
			pointNeighbor[b2] = pointNeighbor_b2;

			//update pointNeighbor b3
			vector<int> pointNeighbor_b3;
			for (int i = 0; i < pointNeighbor[b3].size() / 2; i++) {
				int b31 = pointNeighbor[b3][2 * i];
				int b32 = pointNeighbor[b3][2 * i + 1];
				if (b31 == b1 && b32 == b2) {
					pointNeighbor_b3.push_back(b1);
					pointNeighbor_b3.push_back(index);
					pointNeighbor_b3.push_back(index);
					pointNeighbor_b3.push_back(b2);

					pointNeighbor[index].push_back(b3);
					pointNeighbor[index].push_back(b1);
					pointNeighbor[index].push_back(b2);
					pointNeighbor[index].push_back(b3);

				}
				else if (b31 == b2 && b32 == b1) {
					pointNeighbor_b3.push_back(b2);
					pointNeighbor_b3.push_back(index);
					pointNeighbor_b3.push_back(index);
					pointNeighbor_b3.push_back(b1);

					pointNeighbor[index].push_back(b3);
					pointNeighbor[index].push_back(b2);
					pointNeighbor[index].push_back(b1);
					pointNeighbor[index].push_back(b3);

				}
				else {
					pointNeighbor_b3.push_back(b31);
					pointNeighbor_b3.push_back(b32);
				}
			}
			pointNeighbor[b3].clear();
			pointNeighbor[b3] = pointNeighbor_b3;
			return true;
		}
		else {
			return false;
		}
	}

	vector<Point_In_Edge> TV_check_Edge(int b1, int b2) {

		if (b1 == 11660 && b2 == 24072) {

			cout << "Hello!" << endl;

		}

		vector<int> storePoint;
		vector<int> detectPoint;
		for (int i = 0; i < pointEdge[b1].size(); i++) {
			int b1n = pointEdge[b1][i];
			if (b1n == b2) {
				continue;
			}
			if (!TC_exist_Edge(b1n, storePoint)) {
				storePoint.push_back(b1n);
				if (TC_check_SubEdge_In_Edge(b1, b2, b1n)) {
					detectPoint.push_back(b1n);
				}
				else {
					continue;
				}
			}
		}
		for (int i = 0; i < pointEdge[b2].size(); i++) {
			int b2n = pointEdge[b2][i];
			if (b2n == b1) {
				continue;
			}
			if (!TC_exist_Edge(b2n, storePoint)) {
				storePoint.push_back(b2n);
				if (TC_check_SubEdge_In_Edge(b1, b2, b2n)) {
					detectPoint.push_back(b2n);
				}
				else {
					continue;
				}
			}
		}
		vector<Point_In_Edge> pieList;
		for (int i = 0; i < detectPoint.size(); i++) {

			Point_In_Edge pie;
			pie.b1 = b1;
			pie.b2 = b2;
			pie.index = detectPoint[i];
			pieList.push_back(pie);

		}
		return pieList;

	}

	bool TC_check_SubEdge_In_Edge(int b1, int b2, int checkIndex) {

		if (TC_TriangularProcessed(b1, b2, checkIndex) != -1) {
			return false;
		}

		vector<double> p1 = pointSet[b1];
		vector<double> p2 = pointSet[b2];
		vector<double> checkP = pointSet[checkIndex];

		//judge the linear relationship
		vector<double> v1(3);
		v1[0] = checkP[0] - p1[0];
		v1[1] = checkP[1] - p1[1];
		v1[2] = checkP[2] - p1[2];

		vector<double> v2(3);
		v2[0] = checkP[0] - p2[0];
		v2[1] = checkP[1] - p2[1];
		v2[2] = checkP[2] - p2[2];

		vector<double> v12(3);
		v12[0] = p2[0] - p1[0];
		v12[1] = p2[1] - p1[1];
		v12[2] = p2[2] - p1[2];

		vector<double> v21(3);
		v21[0] = -v12[0];
		v21[1] = -v12[1];
		v21[2] = -v12[2];

		double angle1 = TC_VectorAngel(v1, v12);
		double angle2 = TC_VectorAngel(v2, v21);

		if (angle1 < 0.0005 && angle2 < 0.0005) {
			return true;
		}
		else {
			return false;
		}
	}

	double TC_VectorAngel(vector<double> v1, vector<double> v2) {

		double v12 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
		double v1L = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
		double v2L = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
		double cosV12 = v12 / (v1L * v2L);
		if (cosV12 > 1) {
			cosV12 = 1;
		}
		if (cosV12 < -1) {
			cosV12 = -1;
		}
		return acos(cosV12);

	}

	void TC_init_Edge() {

		pointEdge.resize(pointNeighbor.size());
		for (int i = 0; i < pointNeighbor.size(); i++) {
			vector<int> pointEdge_i;
			for (int j = 0; j < pointNeighbor[i].size(); j++) {
				int b2 = pointNeighbor[i][j];
				if (TC_exist_Edge(b2, pointEdge_i)) {
					continue;
				}
				else {
					pointEdge_i.push_back(b2);
				}
			}
			pointEdge[i] = pointEdge_i;
		}

	}

	bool TC_exist_Edge(int b1, vector<int> pointEdge_i) {

		for (int i = 0; i < pointEdge_i.size(); i++) {
			if (b1 == pointEdge_i[i]) {
				return true;
			}
		}
		return false;
	}

	int TC_TriangularProcessed(int p1, int p2, int p3) {

		for (int i = 0; i < pointNeighbor[p1].size() / 2; i++) {

			int b2 = pointNeighbor[p1][2 * i];
			int b3 = pointNeighbor[p1][2 * i + 1];
			if (b2 == p2 && b3 == p3) {
				return i;
			}
			else if (b2 == p3 && b3 == p2) {
				return i;
			}
			else {
				continue;
			}
		}

		//cout << "Error! the triangular can not be found!" << endl;
		return -1;

	}

	void TC_Print() {

		cout << "TC_Print:" << endl;

		cout << "2320: ";
		for (int i = 0; i < pointNeighbor[2320].size(); i++) {
			if (i == pointNeighbor[2320].size() - 1) {
				cout << pointNeighbor[2320][i] << ".";
			}
			else {
				cout << pointNeighbor[2320][i] << ",";
			}
		}

		cout << endl;

		cout << "2319: ";
		for (int i = 0; i < pointNeighbor[2319].size(); i++) {
			if (i == pointNeighbor[2319].size() - 1) {
				cout << pointNeighbor[2319][i] << ".";
			}
			else {
				cout << pointNeighbor[2319][i] << ",";
			}
		}

		cout << endl;

		cout << "2318: ";
		for (int i = 0; i < pointNeighbor[2318].size(); i++) {
			if (i == pointNeighbor[2318].size() - 1) {
				cout << pointNeighbor[2318][i] << ".";
			}
			else {
				cout << pointNeighbor[2318][i] << ",";
			}
		}

		cout << endl;


	}

	//final round
	void TC_Neighbor2Face() {

		faceSet.clear();

		vector<vector<bool>> pointNeighbor_Judge(pointNeighbor.size());
		for (int i = 0; i < pointNeighbor.size(); i++) {
			vector<bool> pointNeighbor_Judge_i(pointNeighbor[i].size() / 2, false);
			pointNeighbor_Judge[i] = pointNeighbor_Judge_i;
		}

		for (int i = 0; i < pointNeighbor_Judge.size(); i++) {
			int b1 = i;
			vector<int> faceSetNew_i(3);
			faceSetNew_i[0] = b1;
			for (int j = 0; j < pointNeighbor_Judge[i].size(); j++) {
				if (!pointNeighbor_Judge[i][j]) {
					int b2 = pointNeighbor[i][2 * j];
					int b3 = pointNeighbor[i][2 * j + 1];

					if (b1 == b2 || b1 == b3 || b2 == b3) {
						cout << "error edge!" << endl;
					}

					faceSetNew_i[1] = b2;
					faceSetNew_i[2] = b3;
					faceSet.push_back(faceSetNew_i);
					int b2Num = TC_TriangularProcessed(b2, b3, b1);//set have processed triangular to true
					if (b2Num == -1) {
						cout << "Hello!";
					}
					pointNeighbor_Judge[b2][b2Num] = true;
					int b3Num = TC_TriangularProcessed(b3, b1, b2);
					if (b3Num == -1) {
						cout << "Hello!";
					}
					pointNeighbor_Judge[b3][b3Num] = true;
				}
			}
		}

	}

	void TC_SaveOBJ(string fileName, vector<vector<double>> points, vector<vector<int>> facet) {

		ofstream f1(fileName);

		//f1 << points.size() << " " << facet.size() << " " << 0 << endl;

		for (int i = 0; i < points.size(); i++) {
			f1 << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2] << " " << endl;
		}

		for (int i = 0; i < facet.size(); i++) {
			f1 << "f " << facet[i][0] + 1 << " " << facet[i][1] + 1 << " " << facet[i][2] + 1 << endl;
		}

		f1.close();

	}



};
