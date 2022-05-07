/*********************************************************************************

			  Mesh Optimization to achieve isotropic remeshing (Cuda)

						    Updating in 2021/11/02

						      By Dr. Chenlei Lv

		The functions includes:
		Cropping mesh
		[Botsch 2004] A Remeshing Approach to Multiresolution Modeling.

*********************************************************************************/
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h> 
#include <algorithm>
#include <pcl/io/pcd_io.h>  //File input$output
#include <pcl/octree/octree_search.h>  //octree define
#include <pcl/point_types.h>  //point type
#include <pcl/kdtree/kdtree_flann.h>
#include "Mesh_Geometric.cpp"
using namespace std;

extern "C" void AnIsotropic_Unit_Cuda(

	std::vector<vector<double>> points,//source point cloud with different centers	
	std::vector<vector<int>> pointNeighbor,//parameter of target point cloud voxel structure: minXYZ,maxXYZ,unitSize
	vector<double> cu_ave,//parameter of target point cloud voxel structure: XYZNumber
	vector<int> points_Keep,
	double** pointsNew,
	int** pointNeighborNew

);

class AnIsotropic_Cuda {

private:

	int originalNumber;
	string fileName;
	vector<double> cu_ave;
	vector<vector<double>> points;//point cloud
	vector<vector<int>> faceInfor;//trangulars
	vector<vector<int>> pointNeighbor;//neighbor faces
	vector<int> points_Keep;//some points should not be moved or changed
	double L_ave; //the ave edge length of all borders
	vector<double> L_range;	
	int optIter;

public:

	void AnIsotropic_Cuda_init(string fileN, int adtParameter, vector<vector<double>> pointsi, vector<vector<int>> faceInfori) {

		originalNumber = pointsi.size();
		fileName = fileN;
		cout << "Mesh optimization start!" << endl;
		points = pointsi;
		faceInfor = faceInfori;//the index of first point is 1		
		//init neighbor structure, correct point index to first point index is 0
		pointNeighbor.resize(points.size());

		for (int i = 0; i < faceInfor.size(); i++) {

			int b1 = faceInfor[i][0];
			int b2 = faceInfor[i][1];
			int b3 = faceInfor[i][2];
			pointNeighbor[b1].push_back(b2);
			pointNeighbor[b1].push_back(b3);
			pointNeighbor[b2].push_back(b3);
			pointNeighbor[b2].push_back(b1);
			pointNeighbor[b3].push_back(b1);
			pointNeighbor[b3].push_back(b2);

		}

		//init ave length of trangular border
		MeshOptimization_AveEdge_init();
		cout << "L_ave:" << L_ave << endl;

		cout << "Compute Curvature:" << endl;
		MeshOptimization_Apt_L_init(adtParameter);		

		MeshOptimization_EdgePointCheck();
		cout << "Edge Point:" << points_Keep.size() << endl;		

	}

	void AnIsotropic_Cuda_Start(double multiple, int iter, bool judge) {

		L_ave = L_ave * multiple;

		for (int i = 0; i < cu_ave.size(); i++) {
			cu_ave[i] = cu_ave[i] * multiple;		
		}

		optIter = iter;				
		double** pointsNew;
		int** pointNeighborNew;

		AnIsotropic_Unit_Cuda(
			points,//source point cloud with different centers	
			pointNeighbor,//parameter of target point cloud voxel structure: minXYZ,maxXYZ,unitSize
			cu_ave,//parameter of target point cloud voxel structure: XYZNumber
			points_Keep,
			pointsNew,
			pointNeighborNew
		);

		MeshOptimization_UpdateStructure_Face();		
		cout << endl;

	}

	vector<vector<double>> MeshOptimization_GetPoints() {
		return points;
	}

	vector<vector<int>> MeshOptimization_GetFaces() {
		return faceInfor;
	}

private:
	
	void MeshOptimization_Apt_L_init(int classIndex) {

		cu_ave.clear();
		MeshGeometric mg;
		mg.MeshGeometric_init_Harmonic(points, pointNeighbor);
		cu_ave = mg.MeshGeometric_Get_HN();
				
		//mg.MeshGeometric_init(points, pointNeighbor);
		//cu_ave = mg.MeshGeometric_Get_N();

		//average
		int step = 3;
		while (step--) {
			vector<double> pL_ave_store = cu_ave;
			for (int i = 0; i < points.size(); i++) {
				vector<int> pointsN = pointNeighbor[i];
				vector<int> pointsNRerepeat;
				vector<double> weight_list;
				double weight_Sum = 0;
				for (int j = 0; j < pointsN.size(); j++) {
					int indexij = pointsN[j];
					if (MeshOptimization_Exist(indexij, pointsNRerepeat)) {
						continue;
					}
					else {
						pointsNRerepeat.push_back(indexij);
						double weightij = MeshOptimization_EdgeLength(i, indexij);
						weight_Sum = weight_Sum + weightij;
						weight_list.push_back(weightij);
					}
				}
				double L_ave_i = 0;
				for (int j = 0; j < pointsNRerepeat.size(); j++) {
					L_ave_i = L_ave_i + (weight_list[j] / weight_Sum) * pL_ave_store[pointsNRerepeat[j]];
				}
				cu_ave[i] = L_ave_i;
			}
		}

		double disMax = -9999;
		double disMin = 9999;
		for (int i = 0; i < cu_ave.size(); i++) {
			if (cu_ave[i] < disMin) {
				disMin = cu_ave[i];
			}
			if (cu_ave[i] > disMax) {
				disMax = cu_ave[i];
			}
		}
		vector<double> cu_ave_t = cu_ave;
		std::sort(cu_ave_t.begin(), cu_ave_t.end());
		//different adaptively isotropic remeshing method
		
		if (classIndex == 1) {
			double li1 = cu_ave_t[cu_ave_t.size() / 5];
			double li2 = cu_ave_t[cu_ave_t.size() * 2 / 5];
			double li3 = cu_ave_t[cu_ave_t.size() * 3 / 5];
			double li4 = cu_ave_t[cu_ave_t.size() * 4 / 5];
			for (int i = 0; i < cu_ave.size(); i++) {
				if (cu_ave[i] < li1) {
					cu_ave[i] = 2 * L_ave;
				}
				else if (cu_ave[i] < li2 && cu_ave[i] >= li1) {
					cu_ave[i] = 1.5 * L_ave;
				}
				else if (cu_ave[i] < li3 && cu_ave[i] >= li2) {
					cu_ave[i] = 1.0 * L_ave;
				}
				else if (cu_ave[i] < li4 && cu_ave[i] >= li3) {
					cu_ave[i] = 0.8 * L_ave;
				}
				else {
					cu_ave[i] = 0.6 * L_ave;
				}
			}
		}
		else if (classIndex == 2) {
			double li1 = cu_ave_t[cu_ave_t.size() / 8];
			double li2 = cu_ave_t[cu_ave_t.size() * 2 / 8];
			double li3 = cu_ave_t[cu_ave_t.size() * 3 / 8];
			double li4 = cu_ave_t[cu_ave_t.size() * 4 / 8];
			double li5 = cu_ave_t[cu_ave_t.size() * 5 / 8];
			double li6 = cu_ave_t[cu_ave_t.size() * 6 / 8];
			double li7 = cu_ave_t[cu_ave_t.size() * 7 / 8];

			for (int i = 0; i < cu_ave.size(); i++) {
				if (cu_ave[i] < li1) {
					cu_ave[i] = 3 * L_ave;
				}
				else if (cu_ave[i] < li2 && cu_ave[i] >= li1) {
					cu_ave[i] = 2.5 * L_ave;
				}
				else if (cu_ave[i] < li3 && cu_ave[i] >= li2) {
					cu_ave[i] = 2 * L_ave;
				}
				else if (cu_ave[i] < li4 && cu_ave[i] >= li3) {
					cu_ave[i] = 1.5 * L_ave;
				}
				else if (cu_ave[i] < li5 && cu_ave[i] >= li4) {
					cu_ave[i] = 1.2 * L_ave;
				}
				else if (cu_ave[i] < li6 && cu_ave[i] >= li5) {
					cu_ave[i] = 1 * L_ave;
				}
				else if (cu_ave[i] < li7 && cu_ave[i] >= li6) {
					cu_ave[i] = 0.8 * L_ave;
				}
				else {
					cu_ave[i] = 0.7 * L_ave;
				}
			}
		}
		else {

		}
	}
	
	//Finally output the 
	void MeshOptimization_UpdateStructure_Face() {

		vector<vector<int>> pointNeighborNew = pointNeighbor;
		vector<vector<int>> faceInforNew;
		//update face list
		for (int i = 0; i < pointNeighborNew.size(); i++) {
			int b1 = i;
			vector<int> b1_N = pointNeighborNew[b1];
			for (int j = 0; j < pointNeighborNew[i].size() / 2; j++) {
				int b2 = pointNeighborNew[i][2 * j];
				int b3 = pointNeighborNew[i][2 * j + 1];
				if (b2 == -1 || b3 == -1) {
					continue;
				}
				else {
					//f1 << "f " << i + 1 << " " << pointNeighborNew[i][2 * j] + 1 << " " << pointNeighborNew[i][2 * j + 1] + 1 << endl;
					vector<int> faceInforNew_i;
					faceInforNew_i.push_back(i);
					faceInforNew_i.push_back(pointNeighborNew[i][2 * j]);
					faceInforNew_i.push_back(pointNeighborNew[i][2 * j + 1]);
					faceInforNew.push_back(faceInforNew_i);
					for (int k = 0; k < pointNeighborNew[b2].size() / 2; k++) {
						int b2n2 = pointNeighborNew[b2][2 * k];
						int b2n3 = pointNeighborNew[b2][2 * k + 1];
						if ((b2n2 == b1 || b2n2 == b3) && (b2n3 == b1 || b2n3 == b3)) {
							pointNeighborNew[b2][2 * k] = -1;
							pointNeighborNew[b2][2 * k + 1] = -1;
						}
					}
					for (int k = 0; k < pointNeighborNew[b3].size() / 2; k++) {
						int b3n2 = pointNeighborNew[b3][2 * k];
						int b3n3 = pointNeighborNew[b3][2 * k + 1];
						if ((b3n2 == b1 || b3n2 == b2) && (b3n3 == b1 || b3n3 == b2)) {
							pointNeighborNew[b3][2 * k] = -1;
							pointNeighborNew[b3][2 * k + 1] = -1;
						}
					}
				}
			}
		}
		faceInfor.clear();
		faceInfor = faceInforNew;
	}

	//*******************************************Basic Function*****************************************************
	//save mesh
	void MeshOptimization_SaveOBJ() {

		string fileName = "Remesh\\optimization\\" + fileName + ".obj";

		ofstream f1(fileName);
		//ofstream f1(fileName, ios::app);//write in the end

		//f1 << points.size() << " " << facet.size() << " " << 0 << endl;

		for (int i = 0; i < points.size(); i++) {
			f1 << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2] << " " << endl;
		}

		for (int i = 0; i < pointNeighbor.size(); i++) {
			int b1 = i;
			vector<int> b1_N = pointNeighbor[b1];
			for (int j = 0; j < pointNeighbor[i].size() / 2; j++) {
				int b2 = pointNeighbor[i][2 * j];
				int b3 = pointNeighbor[i][2 * j + 1];
				if (b2 == -1 || b3 == -1) {
					continue;
				}
				else {
					f1 << "f " << i + 1 << " " << pointNeighbor[i][2 * j] + 1 << " " << pointNeighbor[i][2 * j + 1] + 1 << endl;
					for (int k = 0; k < pointNeighbor[b2].size() / 2; k++) {
						int b2n2 = pointNeighbor[b2][2 * k];
						int b2n3 = pointNeighbor[b2][2 * k + 1];
						if ((b2n2 == b1 || b2n2 == b3) && (b2n3 == b1 || b2n3 == b3)) {
							pointNeighbor[b2][2 * k] = -1;
							pointNeighbor[b2][2 * k + 1] = -1;
						}
					}
					for (int k = 0; k < pointNeighbor[b3].size() / 2; k++) {
						int b3n2 = pointNeighbor[b3][2 * k];
						int b3n3 = pointNeighbor[b3][2 * k + 1];
						if ((b3n2 == b1 || b3n2 == b2) && (b3n3 == b1 || b3n3 == b2)) {
							pointNeighbor[b3][2 * k] = -1;
							pointNeighbor[b3][2 * k + 1] = -1;
						}
					}
				}
			}
		}

		f1.close();


	}

	void MeshOptimization_SaveOBJ_Face() {

		string fin = "Data\\subRemesh\\" + fileName + "_AdRemesh.obj";

		ofstream f1(fin);

		//f1 << points.size() << " " << facet.size() << " " << 0 << endl;

		for (int i = 0; i < points.size(); i++) {
			f1 << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2] << " " << endl;
		}

		for (int i = 0; i < faceInfor.size(); i++) {
			f1 << "f " << faceInfor[i][0] + 1 << " " << faceInfor[i][1] + 1 << " " << faceInfor[i][2] + 1 << endl;
		}

		f1.close();


	}

	//compute ave edge length of a mesh
	void MeshOptimization_AveEdge_init() {//computer the ave edge l for trangular edge estimation

		//select 1000 trangulars for estimation

		double l_sum = 0;
		for (int i = 0; i < faceInfor.size(); i++) {
			int b1 = faceInfor[i][0];
			int b2 = faceInfor[i][1];
			int b3 = faceInfor[i][2];
			double l1 = MeshOptimization_EdgeLength(b1, b2);
			double l2 = MeshOptimization_EdgeLength(b2, b3);
			double l3 = MeshOptimization_EdgeLength(b3, b1);
			double l_ave_i = (l1 + l2 + l3) / 3;
			l_sum = l_sum + l_ave_i;
		}
		L_ave = l_sum / faceInfor.size();
	}

	//return the length of an edge (b1, b2)
	double MeshOptimization_EdgeLength(int b1, int b2) {

		vector<double> p1 = points[b1];
		vector<double> p2 = points[b2];

		double lengthE = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
			(p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]));

		return lengthE;

	}

	//judge the point is in a list
	bool MeshOptimization_Exist(int b, vector<int> v) {

		for (int i = 0; i < v.size(); i++) {
			if (b == v[i]) {
				return true;
			}
		}
		return false;
	}

	int MeshOptimization_Exist_Index(int b, vector<int> v) {

		for (int i = 0; i < v.size(); i++) {
			if (b == v[i]) {
				return i;
			}
		}
		return -1;
	}
	
	//return the Valence value of a edge (p1, p2). 1: edge; 2: normal; 3: invalid
	int MeshOptimization_Edge_Valence(int p1, int p2) {
		int valence = 0;
		for (int i = 0; i < pointNeighbor[p1].size() / 2; i++) {
			if (pointNeighbor[p1][2 * i] == p2 || pointNeighbor[p1][2 * i + 1] == p2) {
				valence++;
			}
		}
		return valence;
	}

	//return if a point is an edge point
	bool MeshOptimization_Edge_Point(int p1) {

		for (int i = 0; i < pointNeighbor[p1].size(); i++) {
			int p2 = pointNeighbor[p1][i];
			int valenceE = MeshOptimization_Edge_Valence(p1, p2);
			if (valenceE == 1) {
				return true;
			}
		}
		return false;
	}
	
	void MeshOptimization_EdgePointCheck() {

		points_Keep.clear();
		for (int i = 0; i < points.size(); i++) {
			if (MeshOptimization_Edge_Point(i)) {
				//pointNew[i] = pi;
				points_Keep.push_back(i);
			}
		}

	}
	
};


