/*********************************************************************************

							  Sub-Remeshing Scheme

							 Updating in 2021/05/19

							   By Dr. Chenlei Lv

			The functions includes:
			1. Load Mesh;
			2. Sub-remeshing scheme processing;
			3. Some basic vector computation

**********************************************************************************/
#pragma once
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
#pragma once


class Point_In {
public:
	int b1;
	int b2;
	int index;
	vector<double> position;
};

class Sub_Remeshing {

private:

	vector<vector<double>> pointSet;
	vector<vector<int>> faceSet;
	vector<vector<int>> faceSetNew; // new face infor
	vector<vector<int>> pointNeighbor;//store the triangular face 
	vector<vector<int>> pointEdge;//store the edge	

	//new point insert into the V_list
	vector<vector<Point_In>> pointInsert;

	//update point and face information	
	vector<vector<int>> faceNewInsert;

	//point number after edge split
	int indexGlobal;

	//limit edge length
	double E_ave;

public:

	void Sub_Remeshing_init(vector<vector<double>> pointSet_input, vector<vector<int>> faceSet_input) {

		cout << "Sub_Remeshing init run!" << endl;
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
		//Remove repeat item in neighbor structure
		Sub_Removerepeat_Neighbor();

		//init edge vector
		Sub_init_Edge();

		//compute average length;
		Sub_Compute_E_ave();
	}

	void Sub_Remeshing_init(vector<vector<double>> pointSet_input, vector<vector<int>> faceSet_input, double multi) {

		cout << "Sub_Remeshing init run!" << endl;
		pointSet = pointSet_input;
		faceSet = faceSet_input;
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
		//Remove repeat item in neighbor structure
		Sub_Removerepeat_Neighbor();

		//init edge vector
		Sub_init_Edge();

		//compute average length;
		Sub_Compute_E_ave();

		E_ave = E_ave * multi;
	}

	void Sub_Remeshing_Start(string fileNameModel) {

		cout << "Sub_Remeshing Start:" << endl;
		vector<vector<int>> E_distorted = Sub_SearchEdge();
		Sub_ReconnectFace(E_distorted);
		Sub_FrashData();

		cout << "New Point:" << pointSet.size() << endl;
		cout << "New Face:" << faceSetNew.size() << endl;

		//string obj_Sub_Remesh = "Data\\subRemesh\\" + fileNameModel + "_sub.obj";
		//Sub_SaveOBJ(obj_Sub_Remesh, pointSet, faceSetNew);

		faceSet.clear();
		faceSet = faceSetNew;

		cout << "Sub_Remeshing finish." << endl;

	}

	vector<vector<double>> Sub_Get_PointSet() {
		return pointSet;
	}

	vector<vector<int>> Sub_Get_FaceSet() {
		return faceSet;
	}

private:

	//*************************************init*****************************************
	void Sub_Removerepeat_Neighbor() {

		for (int i = 0; i < pointNeighbor.size(); i++) {

			vector<int> pointNeighbor_i;
			for (int j = 0; j < pointNeighbor[i].size() / 2; j++) {
				int b2 = pointNeighbor[i][2 * j];
				int b3 = pointNeighbor[i][2 * j + 1];
				if (Sub_exist_Neighbor(b2, b3, pointNeighbor_i)) {
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

	void Sub_init_Edge() {

		pointEdge.resize(pointNeighbor.size());
		for (int i = 0; i < pointNeighbor.size(); i++) {
			vector<int> pointEdge_i;
			for (int j = 0; j < pointNeighbor[i].size(); j++) {
				int b2 = pointNeighbor[i][j];
				if (Sub_exist_Edge(b2, pointEdge_i)) {
					continue;
				}
				else {
					pointEdge_i.push_back(b2);
				}
			}
			pointEdge[i] = pointEdge_i;
		}

	}

	vector<vector<int>> Sub_SearchEdge() {

		vector<vector<int>> E_distorted(pointEdge.size());
		pointInsert.resize(pointEdge.size());

		int indexTem = pointEdge.size();

		int recentIndex = pointEdge.size();
		for (int i = 0; i < pointEdge.size(); i++) {
			int b1 = i;
			vector<Point_In> Pii;
			vector<int> E_distorted_i;
			for (int j = 0; j < pointEdge[i].size(); j++) {
				int b2 = pointEdge[i][j];
				double dis_ij = Sub_Distance(b1, b2);
				if (dis_ij > 2 * E_ave && b1 < b2) {
					E_distorted_i.push_back(b2);
					indexTem = Sub_Insert_Point_In(b1, b2, indexTem, dis_ij);
				}
			}
			E_distorted[i] = E_distorted_i;
		}

		indexGlobal = indexTem;

		return E_distorted;

	}

	//insert new points into the pointInsert, and return the point index for next round insert. 
	int Sub_Insert_Point_In(int b1, int b2, int index, double distance12) {

		int indexTem = index;
		vector<double> pb1 = pointSet[b1];
		vector<double> pb2 = pointSet[b2];
		int insertNum = distance12 / E_ave;
		int insertNumF = insertNum + 1;
		for (int i = 0; i < insertNum; i++) {
			Point_In pi;
			pi.b1 = b1;
			pi.b2 = b2;
			pi.index = indexTem;
			indexTem++;
			double rate = (double)(i + 1) / (double)insertNumF;
			vector<double> piPosition(3);
			piPosition[0] = pb1[0] * (1 - rate) + pb2[0] * rate;
			piPosition[1] = pb1[1] * (1 - rate) + pb2[1] * rate;
			piPosition[2] = pb1[2] * (1 - rate) + pb2[2] * rate;
			pi.position = piPosition;
			pointInsert[b1].push_back(pi);
		}
		return indexTem;

	}

	//reconnect face
	void Sub_ReconnectFace(vector<vector<int>> E_distorted) {

		//define trigular judge vector
		vector<vector<bool>> pointNeighbor_Judge(pointNeighbor.size());
		for (int i = 0; i < pointNeighbor.size(); i++) {

			vector<bool> pointNeighbor_Judge_i(pointNeighbor[i].size() / 2, false);
			pointNeighbor_Judge[i] = pointNeighbor_Judge_i;

		}

		//start to remesh the trigular faces		
		int count = 0;
		for (int i = 0; i < pointNeighbor.size(); i++) {
			if (i % 100 == 0) {
				cout << pointNeighbor.size() - i << ";";
			}
			int b1 = i;
			for (int j = 0; j < pointNeighbor[i].size() / 2; j++) {
				int b2 = pointNeighbor[i][2 * j];
				int b3 = pointNeighbor[i][2 * j + 1];
				if (pointNeighbor_Judge[i][j]) {
					continue;
				}
				else {
					//cout << "count:" << count << endl;
					//count++;
					//if (count == 363) {
						//cout << "Hello" << endl;					
					//}
					//extract edge, b1, b2
					vector<vector<int>> faceNew = Sub_Insert_Face(b1, b2, b3);
					if (faceNew.size() != 0) {
						faceNewInsert.insert(faceNewInsert.end(), faceNew.begin(), faceNew.end());
						//set the same faces stored in other points to true
						int b2Num = Sub_TriangularProcessed(b2, b3, b1);//set have processed triangular to true
						pointNeighbor_Judge[b2][b2Num] = true;
						int b3Num = Sub_TriangularProcessed(b3, b1, b2);
						pointNeighbor_Judge[b3][b3Num] = true;
						pointNeighbor_Judge[i][j] = true;
					}
				}
			}
		}

		cout << endl;

		for (int i = 0; i < pointNeighbor_Judge.size(); i++) {
			int b1 = i;
			vector<int> faceSetNew_i(3);
			faceSetNew_i[0] = b1;
			for (int j = 0; j < pointNeighbor_Judge[i].size(); j++) {
				if (!pointNeighbor_Judge[i][j]) {
					int b2 = pointNeighbor[i][2 * j];
					int b3 = pointNeighbor[i][2 * j + 1];
					faceSetNew_i[1] = b2;
					faceSetNew_i[2] = b3;
					faceSetNew.push_back(faceSetNew_i);
					int b2Num = Sub_TriangularProcessed(b2, b3, b1);//set have processed triangular to true
					pointNeighbor_Judge[b2][b2Num] = true;
					int b3Num = Sub_TriangularProcessed(b3, b1, b2);
					pointNeighbor_Judge[b3][b3Num] = true;
				}
			}
		}
		cout << "faceSetNew:" << faceSetNew.size() << endl;
	}

	//return new face list
	vector<vector<int>> Sub_Insert_Face(int b1, int b2, int b3) {

		//check if require insert point
		vector<vector<int>> faceVector;
		if (pointInsert[b1].size() == 0 && pointInsert[b2].size() == 0 && pointInsert[b3].size() == 0) {
			return faceVector;
		}

		//store the original face infor
		vector<int> faceOriginal;
		faceOriginal.push_back(b1);
		faceOriginal.push_back(b2);
		faceOriginal.push_back(b3);

		//chech three edge if they take new insert point
		//achieve longest edge:
		double d12 = Sub_Distance(b1, b2);
		double d23 = Sub_Distance(b2, b3);
		double d13 = Sub_Distance(b1, b3);

		int startPoint;
		int endPoint;
		int middlePoint;

		if (d12 >= d23 && d12 >= d13) {
			middlePoint = b3;
			if (d23 > d13) {
				startPoint = b2;
				endPoint = b1;
			}
			else {
				startPoint = b1;
				endPoint = b2;
			}
		}
		else if (d23 >= d12 && d23 >= d13) {
			middlePoint = b1;
			if (d12 > d13) {
				startPoint = b2;
				endPoint = b3;
			}
			else {
				startPoint = b3;
				endPoint = b2;
			}
		}
		else {
			middlePoint = b2;
			if (d12 > d23) {
				startPoint = b1;
				endPoint = b3;
			}
			else {
				startPoint = b3;
				endPoint = b1;
			}
		}

		vector<Point_In> b1pV; //longest edge;
		vector<Point_In> b2pV; //opposite edge to longest edge;
		vector<Point_In> b3pV; //third edge;

		//insert b1pV
		b1pV = Sub_ExactInsertPoint(startPoint, endPoint);
		b2pV = Sub_ExactInsertPoint(startPoint, middlePoint);
		b3pV = Sub_ExactInsertPoint(middlePoint, endPoint);

		if (b1pV.size() == 0 && b2pV.size() == 0 && b3pV.size() == 0) {
			return faceVector;
		}

		//start reconnect: b1pV, b2pV, b3pV.
		int b1Num = b1pV.size();
		int b2Num = b2pV.size();
		int b3Num = b3pV.size();

		int b2b3Num = b2Num + b3Num + 1;

		double dL = Sub_Distance(startPoint, endPoint);//
		double dC = Sub_Distance(startPoint, middlePoint);
		double dM = Sub_Distance(endPoint, middlePoint);

		//compute three angles
		//if the triangular face is a Obtuse triangle
		double cosStart = (dL * dL + dC * dC - dM * dM) / (2 * dC * dL); //vertex start
		double cosMiddle = (dM * dM + dC * dC - dL * dL) / (2 * dC * dM); //vertex middle
		double cosEnd = (dM * dM + dL * dL - dC * dC) / (2 * dL * dM); //vertex end

		if (cosStart > 1) {
			cosStart = 1;
		}
		else if (cosStart < -1) {
			cosStart = -1;
		}
		if (cosMiddle > 1) {
			cosMiddle = 1;
		}
		else if (cosMiddle < -1) {
			cosMiddle = -1;
		}
		if (cosEnd > 1) {
			cosEnd = 1;
		}
		else if (cosEnd < -1) {
			cosEnd = -1;
		}


		double angStart = acos(cosStart);
		double angMiddle = acos(cosMiddle);
		double angEnd = acos(cosEnd);

		vector<int> pairPoint;

		if (angMiddle > 1.5) {
			//Obtuse triangle, b1Num, b2b3Num			
			if (b2b3Num > b1Num) {
				for (int i = 1; i <= b2b3Num; i++) {
					int b1NumIndex = (double)b1Num * (double)i / (double)b2b3Num;
					if (b1Num == b1NumIndex) {
						b1NumIndex = b1NumIndex - 1;
					}
					pairPoint.push_back(b1NumIndex);
					pairPoint.push_back(i - 1);
				}
			}
			else {
				for (int i = 1; i <= b1Num; i++) {
					int b23NumIndex = (double)b2b3Num * (double)i / (double)b1Num;
					if (b23NumIndex == b2b3Num) {
						b23NumIndex = b23NumIndex - 1;
					}
					pairPoint.push_back(i - 1);
					pairPoint.push_back(b23NumIndex);
				}
			}
		}
		else {
			for (int i = 1; i <= b1Num; i++) {
				int b2NumIndex = (double)b2Num * (double)i / (double)b1Num;
				pairPoint.push_back(i - 1);
				pairPoint.push_back(b2NumIndex);
			}
		}

		//check pairPoint

		int tbefore = pairPoint[0];
		int tafter = pairPoint[1];
		if (pairPoint.size() >= 4) {
			vector<int> pairPointCheck;
			pairPointCheck.push_back(tbefore);
			pairPointCheck.push_back(tafter);
			for (int i = 1; i < pairPoint.size() / 2; i++) {
				if (pairPoint[2 * i] == tbefore && pairPoint[2 * i + 1] == tafter) {
					continue;
				}
				else {
					pairPointCheck.push_back(pairPoint[2 * i]);
					pairPointCheck.push_back(pairPoint[2 * i + 1]);
					tbefore = pairPoint[2 * i];
					tafter = pairPoint[2 * i + 1];
				}
			}
			pairPoint.clear();
			pairPoint = pairPointCheck;
		}

		vector<Point_In> b1pV_new; //new list
		vector<Point_In> b23pV; //new list

		for (int i = 0; i < pairPoint.size() / 2; i++) {
			int t = pairPoint[2 * i + 1];
			int t0 = pairPoint[2 * i];
			b1pV_new.push_back(b1pV[t0]);

			if (t == b2pV.size()) {
				Point_In inputPoint_i;
				inputPoint_i.b1 = -1;
				inputPoint_i.b2 = -1;
				inputPoint_i.index = middlePoint;
				inputPoint_i.position = pointSet[middlePoint];
				b23pV.push_back(inputPoint_i);
			}
			else if (t > b2pV.size()) {
				b23pV.push_back(b3pV[t - b2pV.size() - 1]);
			}
			else {
				b23pV.push_back(b2pV[t]);
			}
		}

		//add new face
		for (int i = 0; i <= b1pV_new.size(); i++) {
			vector<int> face_i(3);
			vector<Point_In> inputPoint;
			if (i == 0) {
				Point_In inputPoint_i;
				inputPoint_i.b1 = -1;
				inputPoint_i.b2 = -1;
				inputPoint_i.index = startPoint;
				inputPoint_i.position = pointSet[startPoint];
				inputPoint.push_back(inputPoint_i);
				inputPoint.push_back(b1pV_new[0]);
				inputPoint.push_back(b23pV[0]);
				vector<Point_In> newFace = Sub_TriangularNormal(faceOriginal, inputPoint);
				face_i[0] = newFace[0].index;
				face_i[1] = newFace[1].index;
				face_i[2] = newFace[2].index;
				faceVector.push_back(face_i);
			}
			else if (i == b1pV_new.size()) {
				Point_In inputPoint_i;
				inputPoint_i.b1 = -1;
				inputPoint_i.b2 = -1;
				inputPoint_i.index = endPoint;
				inputPoint_i.position = pointSet[endPoint];
				inputPoint.push_back(b1pV_new[i - 1]);
				inputPoint.push_back(b23pV[i - 1]);
				inputPoint.push_back(inputPoint_i);
				vector<Point_In> newFace = Sub_TriangularNormal(faceOriginal, inputPoint);
				face_i[0] = newFace[0].index;
				face_i[1] = newFace[1].index;
				face_i[2] = newFace[2].index;
				faceVector.push_back(face_i);
			}
			else {
				Point_In p1 = b1pV_new[i - 1];
				Point_In p2 = b1pV_new[i];
				Point_In p3 = b23pV[i - 1];
				Point_In p4 = b23pV[i];
				vector<Point_In> plist;
				plist.push_back(p1);
				plist.push_back(p2);
				plist.push_back(p3);
				plist.push_back(p4);
				vector<vector<int>> newFaceIndex = Sub_Triangular_Insert_NewFace(faceOriginal, plist);
				for (int j = 0; j < newFaceIndex.size(); j++) {
					faceVector.push_back(newFaceIndex[j]);
				}
			}
		}
		return faceVector;


	}

	vector<vector<int>> Sub_Triangular_Insert_NewFace(vector<int> Regular, vector<Point_In> p) {

		vector<vector<int>> newFace;

		Point_In p1 = p[0];
		Point_In p2 = p[1];
		Point_In p3 = p[2];
		Point_In p4 = p[3];

		if (p1.index == p2.index) {

			vector<Point_In> pList;
			pList.push_back(p1);
			pList.push_back(p3);
			pList.push_back(p4);
			vector<Point_In> pListNew = Sub_TriangularNormal(Regular, pList);
			vector<int> b1;
			b1.push_back(pListNew[0].index);
			b1.push_back(pListNew[1].index);
			b1.push_back(pListNew[2].index);
			newFace.push_back(b1);
			return newFace;

		}
		else if (p3.index == p4.index) {

			vector<Point_In> pList;
			pList.push_back(p1);
			pList.push_back(p2);
			pList.push_back(p3);
			vector<Point_In> pListNew = Sub_TriangularNormal(Regular, pList);
			vector<int> b1;
			b1.push_back(pListNew[0].index);
			b1.push_back(pListNew[1].index);
			b1.push_back(pListNew[2].index);
			newFace.push_back(b1);
			return newFace;

		}
		else {

			double p14 = sqrt((p1.position[0] - p4.position[0]) * (p1.position[0] - p4.position[0]) +
				(p1.position[1] - p4.position[1]) * (p1.position[1] - p4.position[1]) +
				(p1.position[0] - p4.position[0]) * (p1.position[0] - p4.position[0]));
			double p23 = sqrt((p2.position[0] - p3.position[0]) * (p2.position[0] - p3.position[0]) +
				(p2.position[1] - p3.position[1]) * (p2.position[1] - p3.position[1]) +
				(p2.position[0] - p3.position[0]) * (p2.position[0] - p3.position[0]));

			if (p14 < p23) {

				vector<Point_In> pList;
				pList.push_back(p1);
				pList.push_back(p3);
				pList.push_back(p4);
				vector<Point_In> pListNew = Sub_TriangularNormal(Regular, pList);
				vector<int> b1;
				b1.push_back(pListNew[0].index);
				b1.push_back(pListNew[1].index);
				b1.push_back(pListNew[2].index);
				newFace.push_back(b1);

				vector<Point_In> pList2;
				pList2.push_back(p1);
				pList2.push_back(p2);
				pList2.push_back(p4);
				vector<Point_In> pListNew2 = Sub_TriangularNormal(Regular, pList2);
				vector<int> b2;
				b2.push_back(pListNew2[0].index);
				b2.push_back(pListNew2[1].index);
				b2.push_back(pListNew2[2].index);
				newFace.push_back(b2);

				return newFace;
			}
			else {

				vector<Point_In> pList;
				pList.push_back(p1);
				pList.push_back(p2);
				pList.push_back(p3);
				vector<Point_In> pListNew = Sub_TriangularNormal(Regular, pList);
				vector<int> b1;
				b1.push_back(pListNew[0].index);
				b1.push_back(pListNew[1].index);
				b1.push_back(pListNew[2].index);
				newFace.push_back(b1);

				vector<Point_In> pList2;
				pList2.push_back(p2);
				pList2.push_back(p3);
				pList2.push_back(p4);
				vector<Point_In> pListNew2 = Sub_TriangularNormal(Regular, pList2);
				vector<int> b2;
				b2.push_back(pListNew2[0].index);
				b2.push_back(pListNew2[1].index);
				b2.push_back(pListNew2[2].index);
				newFace.push_back(b2);
				return newFace;
			}
		}
	}

	vector<Point_In> Sub_ExactInsertPoint(int b1, int b2) {

		vector<Point_In> result;
		if (b1 < b2) {
			vector<Point_In> p12 = pointInsert[b1];
			for (int i = 0; i < p12.size(); i++) {
				Point_In p12_i = p12[i];
				if (p12_i.b1 == b1 && p12_i.b2 == b2) {
					result.push_back(p12_i);
				}
				else if (p12_i.b2 == b1 && p12_i.b1 == b2) {
					cout << "Sub_ExactInsertPoint: order confuse!";
				}
				else {
					continue;
				}
			}
		}
		else {
			vector<Point_In> resultC;
			vector<Point_In> p21 = pointInsert[b2];
			for (int i = 0; i < p21.size(); i++) {
				Point_In p21_i = p21[i];
				if (p21_i.b1 == b2 && p21_i.b2 == b1) {
					resultC.push_back(p21_i);
				}
				else if (p21_i.b1 == b1 && p21_i.b2 == b2) {
					cout << "Sub_ExactInsertPoint: order confuse!";
				}
				else {
					continue;
				}
			}
			for (int i = resultC.size() - 1; i >= 0; i--) {
				result.push_back(resultC[i]);
			}
			resultC.clear();
		}
		return result;
	}

	vector<Point_In> Sub_TriangularNormal(vector<int> Regular, vector<Point_In> inputPoint) {

		vector<double> Regular_b1 = pointSet[Regular[0]];
		vector<double> Regular_b2 = pointSet[Regular[1]];
		vector<double> Regular_b3 = pointSet[Regular[2]];

		vector<double> v1(3);
		v1[0] = Regular_b2[0] - Regular_b1[0];
		v1[1] = Regular_b2[1] - Regular_b1[1];
		v1[2] = Regular_b2[2] - Regular_b1[2];

		vector<double> v2(3);
		v2[0] = Regular_b3[0] - Regular_b2[0];
		v2[1] = Regular_b3[1] - Regular_b2[1];
		v2[2] = Regular_b3[2] - Regular_b2[2];

		vector<double> normal_regualr = Sub_Cross(v1, v2);

		vector<double> insert_b1 = inputPoint[0].position;
		vector<double> insert_b2 = inputPoint[1].position;
		vector<double> insert_b3 = inputPoint[2].position;

		vector<double> vi1(3);
		vi1[0] = insert_b2[0] - insert_b1[0];
		vi1[1] = insert_b2[1] - insert_b1[1];
		vi1[2] = insert_b2[2] - insert_b1[2];

		vector<double> vi2(3);
		vi2[0] = insert_b3[0] - insert_b2[0];
		vi2[1] = insert_b3[1] - insert_b2[1];
		vi2[2] = insert_b3[2] - insert_b2[2];

		vector<double> normal_new = Sub_Cross(vi1, vi2);

		bool judge = Sub_VectorAngel(normal_regualr, normal_new);
		if (judge) {
			return inputPoint;
		}
		else {
			vector<Point_In> changeInsert;
			changeInsert.push_back(inputPoint[0]);
			changeInsert.push_back(inputPoint[2]);
			changeInsert.push_back(inputPoint[1]);
			return changeInsert;
		}
	}

	int Sub_TriangularProcessed(int p1, int p2, int p3) {

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

		cout << "Error! the triangular can not be found!" << endl;
		return -1;

	}

	//Cross product
	vector<double> Sub_Cross(vector<double> v1, vector<double> v2) {

		vector<double> normal(3);
		normal[0] = v1[1] * v2[2] - v2[1] * v1[2];
		normal[1] = v2[0] * v1[2] - v2[2] * v1[0];
		normal[2] = v1[0] * v2[1] - v2[0] * v1[1];
		return normal;

	}

	//Vector angel
	bool Sub_VectorAngel(vector<double> v1, vector<double> v2) {

		double v12 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
		double v1L = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
		double v2L = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
		double cosV12 = v12 / (v1L * v2L);

		if (cosV12 > 1) {
			cosV12 = 1;
		}
		else if (cosV12 < -1) {
			cosV12 = -1;
		}
		double angle12 = acos(cosV12);

		if (angle12 < 3.1415 / 2.0) {
			return true;
		}
		else {
			return false;
		}
	}

	void Sub_Compute_E_ave() {//computer average edge length E_ave

		double E_sum = 0;
		int count = 0;
		for (int i = 0; i < pointNeighbor.size(); i++) {
			int p1 = i;
			double E_sum_i = 0;
			for (int j = 0; j < pointNeighbor[i].size(); j++) {
				int p2 = pointNeighbor[i][j];
				double distance12 = Sub_Distance(p1, p2);
				E_sum_i = E_sum_i + distance12;
			}
			if (pointNeighbor[i].size() == 0) {
				count++;
			}
			else {
				E_sum_i = E_sum_i / pointNeighbor[i].size();
				E_sum = E_sum + E_sum_i;
			}
		}
		E_sum = E_sum / (double)(pointNeighbor.size() - count);
		E_ave = E_sum / 2;

	}

	bool Sub_exist_Edge(int b1, vector<int> pointEdge_i) {

		for (int i = 0; i < pointEdge_i.size(); i++) {
			if (b1 == pointEdge_i[i]) {
				return true;
			}
		}
		return false;
	}

	bool Sub_exist_Neighbor(int b2, int b3, vector<int> pointNeighbor_i) {

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

	//compute the distance of two points
	double Sub_Distance(int p1, int p2) {

		vector<double> pp1 = pointSet[p1];
		vector<double> pp2 = pointSet[p2];
		double distance12 = sqrt((pp1[0] - pp2[0]) * (pp1[0] - pp2[0]) +
			(pp1[1] - pp2[1]) * (pp1[1] - pp2[1]) + (pp1[2] - pp2[2]) * (pp1[2] - pp2[2]));
		return distance12;

	}

	void Sub_FrashData() {

		//Exist points
		pointSet;
		faceSet;

		vector<Point_In> list;
		for (int i = 0; i < pointInsert.size(); i++) {
			for (int j = 0; j < pointInsert[i].size(); j++) {
				Point_In pij = pointInsert[i][j];
				list.push_back(pij);
			}
		}

		if (list.size() > 0) {
			int startIndex = list[0].index;
			for (int i = 0; i < list.size(); i++) {
				Point_In pi = list[i];
				if (pi.index == startIndex) {
					pointSet.push_back(pi.position);
					startIndex++;
				}
				else {
					cout << "Error! point index wrong!" << endl;
					break;
				}
			}
		}
		faceSetNew.insert(faceSetNew.end(), faceNewInsert.begin(), faceNewInsert.end());

	}

	void Sub_SaveOBJ(string fileName, vector<vector<double>> points, vector<vector<int>> facet) {

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