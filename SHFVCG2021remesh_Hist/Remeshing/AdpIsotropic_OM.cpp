/*********************************************************************************

			  Mesh Optimization to achieve isotropic remeshing (OpenMP)

						    Updating in 2021/11/03

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

class AnIsotropic_OM {

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
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtreeG;
	int optIter;

public:

	void AnIsotropic_OM_init(string fileN, int adtParameter, vector<vector<double>> pointsi, vector<vector<int>> faceInfori) {

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

		//Kd-tree construct
		kdtreeG = MeshOptimization_KdTreeConstruct(pointsi);

		MeshOptimization_EdgePointCheck();
		cout << "Edge Point:" << points_Keep.size() << endl;		

	}

	void AnIsotropic_OM_Start(double multiple, int iter, bool judge) {

		L_ave = L_ave * multiple;

		for (int i = 0; i < cu_ave.size(); i++) {
			cu_ave[i] = cu_ave[i] * multiple;		
		}

		while (iter) {

			cout << "iter:" << iter << endl;
			//split			
			cout << "split:";
			MeshOptimization_OM_Split();			

			//collapse
			cout << "collapse:";
			MeshOptimization_OM_Collapse();

			//fix structure
			AnIsotropic_OM_CheckErrorBorder();

			//Flip
			cout << "flip:";
			MeshOptimization_OM_Flip();

			//Tangent Smoothing
			cout << "tangent smoothing:";
			MeshOptimization_OM_TangentSmoothing();

			iter--;
		}		

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
	
	//construct kd-tree
	pcl::KdTreeFLANN<pcl::PointXYZ>  MeshOptimization_KdTreeConstruct(vector<vector<double>> seedpoints) {

		pcl::KdTreeFLANN<pcl::PointXYZ> kdtreeSeed;
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->width = seedpoints.size();
		cloud->height = 1;
		cloud->points.resize(cloud->width * cloud->height);
		// fills a PointCloud with random data
		for (int i = 0; i < seedpoints.size(); i++)
		{
			pcl::PointXYZ pxyz;
			cloud->points[i].x = seedpoints[i][0];
			cloud->points[i].y = seedpoints[i][1];
			cloud->points[i].z = seedpoints[i][2];

		}
		kdtreeSeed.setInputCloud(cloud);
		return kdtreeSeed;
	}

	double MeshOptimization_KdTreeSearch(int pindex) {

		vector<int> pointIdxNKNSearch;
		vector<float> pointNKNSquaredDistance;
		int kn = 1;
		pcl::PointXYZ searchPoint;
		searchPoint.x = points[pindex][0];
		searchPoint.y = points[pindex][1];
		searchPoint.z = points[pindex][2];
		kdtreeG.nearestKSearch(searchPoint, kn, pointIdxNKNSearch, pointNKNSquaredDistance);
		int L_index = pointIdxNKNSearch[0];
		return cu_ave[L_index];

	}

	double MeshOptimization_KdTreeSearch(vector<double> pointP) {

		vector<int> pointIdxNKNSearch;
		vector<float> pointNKNSquaredDistance;
		int kn = 1;
		pcl::PointXYZ searchPoint;
		searchPoint.x = pointP[0];
		searchPoint.y = pointP[1];
		searchPoint.z = pointP[2];
		kdtreeG.nearestKSearch(searchPoint, kn, pointIdxNKNSearch, pointNKNSquaredDistance);
		int L_index = pointIdxNKNSearch[0];
		return cu_ave[L_index];

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

	//*******************************************Core Function*****************************************************
//#pragma omp parallel for 

	void MeshOptimization_OM_Split() {

		int countWhile = 10;

		while (1) {

			//1. Check split border
			
			vector<bool> pointJudgeActivate(points.size(), true);//record point can be processed in parallel structure
			vector<vector<int>> Split_list;

			for (int i = 0; i < pointNeighbor.size(); i++) {
				vector<int> b2_list;
				int b1 = i;
				if (!pointJudgeActivate[b1]) {
					continue;				
				}
				for (int j = 0; j < pointNeighbor[i].size(); j++) {
					int b2 = pointNeighbor[i][j];
					if (!pointJudgeActivate[b2]) {
						continue;
					}
					if (!MeshOptimization_Exist(b2, b2_list) && b1 < b2) {
						b2_list.push_back(b2);
						double length12 = MeshOptimization_EdgeLength(b1, b2);
						double L_j1 = MeshOptimization_KdTreeSearch(points[b1]);
						double L_j2 = MeshOptimization_KdTreeSearch(points[b2]);
						if (L_j1 > L_j2) {
							L_j1 = L_j2;
						}
						if (length12 >= 1.33 * L_j1) {
							vector<int> edgeij;
							edgeij.push_back(b1);
							edgeij.push_back(b2);
							Split_list.push_back(edgeij);
							//lock related points
							for (int k = 0; k < pointNeighbor[b1].size(); k++) {
								int b1nk = pointNeighbor[b1][k];
								pointJudgeActivate[b1nk] = false;															
							}
							for (int k = 0; k < pointNeighbor[b2].size(); k++) {
								int b2nk = pointNeighbor[b2][k];
								pointJudgeActivate[b2nk] = false;
							}
						}
					}
				}
			}
			
			cout << Split_list.size();

			//2. Split
			vector<vector<double>> pointNewInsert(Split_list.size());
			vector<vector<int>> pointNewInsertNeighbor(Split_list.size());			
			
			if (countWhile == 0 || Split_list.size() == 0) {
				cout << ".";
				break;
			}
			else {
				cout << ",";
			}
            //parallel processing split    
//#pragma omp parallel for
			for (int i = 0; i < Split_list.size(); i++) {

				int pointsEndIndex = i + points.size();

				int bs1 = Split_list[i][0];
				int bs2 = Split_list[i][1];
				
				//Split: 1 achieve middle point

				vector<double> middlePoint(3);
				middlePoint[0] = (points[bs1][0] + points[bs2][0]) / 2;
				middlePoint[1] = (points[bs1][1] + points[bs2][1]) / 2;
				middlePoint[2] = (points[bs1][2] + points[bs2][2]) / 2;		
				pointNewInsert[i] = middlePoint;

				vector<int> pointNeighbor_b1 = pointNeighbor[bs1];
				vector<int> pointNeighbor_middle;//store the structure of middle point

				//Split: 2 store related points, four neighbot points of the new niddle point.

				vector<int> bs3_list;
				bs3_list.push_back(bs1);
				bs3_list.push_back(bs2);
				for (int i = 0; i < pointNeighbor_b1.size() / 2; i++) {
					int b1i1 = pointNeighbor_b1[2 * i];
					int b1i2 = pointNeighbor_b1[2 * i + 1];
					if (b1i1 == bs2) {
						bs3_list.push_back(b1i2);
					}
					else if (b1i2 == bs2) {
						bs3_list.push_back(b1i1);
					}
					else {
						continue;
					}
				}

				//Split: 3 update the structure	for the related points				

				//vector<vector<int>> T_n;
				for (int ii = 0; ii < bs3_list.size(); ii++) {
					int b_index = bs3_list[ii];
					if (b_index == bs1 || b_index == bs2) {
						for (int j = 0; j < pointNeighbor[b_index].size() / 2; j++) {
							int b_index_j1 = pointNeighbor[b_index][2 * j];
							int b_index_j2 = pointNeighbor[b_index][2 * j + 1];
							if (b_index_j1 == bs1 || b_index_j1 == bs2) {
								pointNeighbor[b_index][2 * j] = pointsEndIndex;
								pointNeighbor_middle.push_back(b_index_j2);
								pointNeighbor_middle.push_back(b_index);
							}
							if (b_index_j2 == bs1 || b_index_j2 == bs2) {
								pointNeighbor[b_index][2 * j + 1] = pointsEndIndex;
								pointNeighbor_middle.push_back(b_index);
								pointNeighbor_middle.push_back(b_index_j1);
							}
						}
					}
					else {
						for (int j = 0; j < pointNeighbor[b_index].size() / 2; j++) {
							int b_index_j1 = pointNeighbor[b_index][2 * j];
							int b_index_j2 = pointNeighbor[b_index][2 * j + 1];
							if ((b_index_j1 == bs1 || b_index_j1 == bs2) && (b_index_j2 == bs1 || b_index_j2 == bs2)) {
								pointNeighbor[b_index][2 * j] = pointsEndIndex;
								pointNeighbor[b_index].push_back(b_index_j1);
								pointNeighbor[b_index].push_back(pointsEndIndex);
								pointNeighbor_middle.push_back(b_index_j2);
								pointNeighbor_middle.push_back(b_index);
								pointNeighbor_middle.push_back(b_index);
								pointNeighbor_middle.push_back(b_index_j1);
							}
						}
					}
				}


				//remove repeat
				vector<int> pointNeighbor_middle_new;
				for (int j = 0; j < pointNeighbor_middle.size() / 2; j++) {

					int a1 = pointNeighbor_middle[2 * j];
					int a2 = pointNeighbor_middle[2 * j + 1];
					if (a1 == -1 || a2 == -1) {
						continue;					
					}
					else {
						pointNeighbor_middle_new.push_back(a1);
						pointNeighbor_middle_new.push_back(a2);
						for (int k = j; k < pointNeighbor_middle.size() / 2; k++) {
							int a1k = pointNeighbor_middle[2 * k];
							int a2k = pointNeighbor_middle[2 * k + 1];
							if (a1 == a1k && a2 == a2k) {
								pointNeighbor_middle[2 * k] = -1;
								pointNeighbor_middle[2 * k + 1] = -1;							
							}
						}					
					}

				}
				pointNewInsertNeighbor[i] = pointNeighbor_middle_new;
			}

			//3. Update point structure
			points.insert(points.end(), pointNewInsert.begin(), pointNewInsert.end());
			pointNeighbor.insert(pointNeighbor.end(), pointNewInsertNeighbor.begin(), pointNewInsertNeighbor.end());

			if (countWhile == 0) {
				break;			
			}
			else {
				countWhile--;			
			}
		
		}

		cout << endl;
	
	}

	void MeshOptimization_OM_Collapse() {

		int countWhile = 10;

		while (1) {

			//1. Check collapse border

			vector<bool> pointJudgeActivate(points.size(), true);//record point can be processed in parallel structure			
			vector<vector<int>> Collapse_List;			
			for (int i = 0; i < pointNeighbor.size(); i++) {					
				vector<int> b2_list;
				int b1 = i;
				if (!pointJudgeActivate[b1]) {
					continue;
				}
				for (int j = 0; j < pointNeighbor[i].size(); j++) {					
					int b2 = pointNeighbor[i][j];
					if (b2 > pointJudgeActivate.size() || b2 < 0) {
						cout << "hello!";					
					}
					if (!pointJudgeActivate[b2]) {
						continue;
					}
					if (!MeshOptimization_Exist(b2, b2_list) && b1 < b2) {
						b2_list.push_back(b2);
						double length12 = MeshOptimization_EdgeLength(b1, b2);
						double L_j1 = MeshOptimization_KdTreeSearch(points[b1]);
						double L_j2 = MeshOptimization_KdTreeSearch(points[b2]);
						if (L_j1 > L_j2) {
							L_j1 = L_j2;
						}
						if (length12 < 0.8 * L_j1) {
							int valence = MeshOptimization_Edge_Valence(b1, b2);
							if (valence <= 2) {//label!
								vector<int> edgeij;
								edgeij.push_back(b1);
								edgeij.push_back(b2);
								if (MeshOptimization_Edge_Point(b1) || MeshOptimization_Edge_Point(b2) || MeshOptimization_Edge_Valence(b1, b2) > 2) {
									continue;
								}
								Collapse_List.push_back(edgeij);
								//lock related points
								for (int k = 0; k < pointNeighbor[b1].size(); k++) {
									int b1nk = pointNeighbor[b1][k];
									pointJudgeActivate[b1nk] = false;
								}
								for (int k = 0; k < pointNeighbor[b2].size(); k++) {
									int b2nk = pointNeighbor[b2][k];
									pointJudgeActivate[b2nk] = false;
								}
							}						
						}
					}
				}
			}

			vector<vector<int>> Collapse_List_new;
			vector<vector<double>> Collapse_List_Middle_new;
			vector<vector<double>> Collapse_List_Middle(Collapse_List.size());
			vector<bool> Collapse_List_Judge(Collapse_List.size(), true);

			//2. Remove collapse border if the new borders are length than regular 
//#pragma omp parallel for
			for (int i = 0; i < Collapse_List.size(); i++) {

				bool Collapse_Invaild = false;
				int b1 = Collapse_List[i][0];
				int b2 = Collapse_List[i][1];
				vector<double> middlePointc(3);
				middlePointc[0] = (points[b1][0] + points[b2][0]) / 2;
				middlePointc[1] = (points[b1][1] + points[b2][1]) / 2;
				middlePointc[2] = (points[b1][2] + points[b2][2]) / 2;
				Collapse_List_Middle[i] = middlePointc;

				if (Collapse_List_Middle[i].size() == 0) {
					cout << "Hello";				
				}

				vector<int> pointNeighbor_b1 = pointNeighbor[b1];
				vector<int> pointNeighbor_b2 = pointNeighbor[b2];
				pointNeighbor_b1.insert(pointNeighbor_b1.end(), pointNeighbor_b2.begin(), pointNeighbor_b2.end());
				for (int j = 0; j < pointNeighbor_b1.size(); j++) {
					int b3 = pointNeighbor_b1[j];
					if (b3 == b1 || b3 == b2) {
						continue;					
					}
					else {						
						double radius12 = sqrt((middlePointc[0] - points[b3][0]) * (middlePointc[0] - points[b3][0]) +
							(middlePointc[1] - points[b3][1]) * (middlePointc[1] - points[b3][1]) +
							(middlePointc[2] - points[b3][2]) * (middlePointc[2] - points[b3][2]));
						double L_j1 = MeshOptimization_KdTreeSearch(middlePointc);
						double L_j2 = MeshOptimization_KdTreeSearch(b3);
						if (L_j1 > L_j2) {
							L_j1 = L_j2;
						}
						if (radius12 > 1.33 * L_j1) {
							Collapse_Invaild = true;
							break;
						}					
					}				
				}
				if (Collapse_Invaild) {
					Collapse_List_Judge[i] = false;
					//break;				
				}
			}

			for (int i = 0; i < Collapse_List.size(); i++) {
				if (Collapse_List_Judge[i]) {
					Collapse_List_new.push_back(Collapse_List[i]);
					Collapse_List_Middle_new.push_back(Collapse_List_Middle[i]);
				}			
			}

			Collapse_List.clear();
			Collapse_List = Collapse_List_new;
			Collapse_List_Middle.clear();
			Collapse_List_Middle = Collapse_List_Middle_new;
			
			cout << Collapse_List.size();
			if (countWhile == 0|| Collapse_List.size() == 0) {
				cout << ".";
				break;
			}
			else {
				cout << ",";
			}
			//3. Collapse; update neighbor structure
			vector<vector<int>> pointNeighbor_middlec(Collapse_List.size());
#pragma omp parallel for
			for (int i = 0; i < Collapse_List.size(); i++) {

				int b1 = Collapse_List[i][0];
				int b2 = Collapse_List[i][1];
				vector<int> pc_list = pointNeighbor[b1];
				pc_list.insert(pc_list.end(), pointNeighbor[b2].begin(), pointNeighbor[b2].end());
				int b_middle_new = points.size() + i;
				for (int j = 0; j < pc_list.size(); j++) {
					if (pc_list[j] == b1 || pc_list[j] == b2|| pc_list[j] == -1) {
						continue;					
					}
					int b3 = pc_list[j];

					for (int k = j; k < pc_list.size(); k++) {

						if (pc_list[k] == b3) {
							pc_list[k] = -1;						
						}					
					}

					for (int k = 0; k < pointNeighbor[b3].size() / 2; k++) {

						int pc1 = pointNeighbor[b3][2 * k];
						int pc2 = pointNeighbor[b3][2 * k + 1];

						if ((pc1 == b1 || pc1 == b2) && (pc2 == b1 || pc2 == b2)) {
							pointNeighbor[b3][2 * k] = -1;
							pointNeighbor[b3][2 * k + 1] = -1;							
						}
						else if (pc1 == b1 || pc1 == b2) {
							pointNeighbor[b3][2 * k] = b_middle_new;
							pointNeighbor_middlec[i].push_back(pc2);
							pointNeighbor_middlec[i].push_back(b3);
						}
						else if (pc2 == b1 || pc2 == b2) {
							pointNeighbor[b3][2 * k + 1] = b_middle_new;
							pointNeighbor_middlec[i].push_back(b3);
							pointNeighbor_middlec[i].push_back(pc1);
						}
						else {
							continue;
						}	
					}
				}			
			}	
		
			//4. Update points and pointNeighbor
			vector<int> pointsIndex(points.size() + Collapse_List_Middle.size());
			for (int i = 0; i < pointsIndex.size(); i++) {
				pointsIndex[i] = i;			
			}
			for (int i = 0; i < Collapse_List.size(); i++) {
				int b1 = Collapse_List[i][0];
				int b2 = Collapse_List[i][1];
				pointsIndex[b1] = -1;
				pointsIndex[b2] = -1;
				for (int j = b1 + 1; j < pointsIndex.size(); j++) {
					if (pointsIndex[j] == -1) {
						continue;					
					}
					else {
						pointsIndex[j] = pointsIndex[j] - 1;					
					}
				}
				for (int j = b2 + 1; j < pointsIndex.size(); j++) {
					if (pointsIndex[j] == -1) {
						continue;
					}
					else {
						pointsIndex[j] = pointsIndex[j] - 1;
					}
				}			
			}

			points.insert(points.end(), Collapse_List_Middle.begin(), Collapse_List_Middle.end());
			pointNeighbor.insert(pointNeighbor.end(), pointNeighbor_middlec.begin(), pointNeighbor_middlec.end());

			for (int i = 0; i < points.size(); i++) {

				if (points[i].size() == 0) {

					cout << "Hello!";
				
				}
			
			}

			vector<vector<double>> points_New;
			vector<vector<int>> pointNeighbor_New;

			//update points and pointNeighbor list;
			for (int i = 0; i < pointsIndex.size(); i++) {

				if (pointsIndex[i] == -1) {
					continue;				
				}
				else {
					vector<double> point_insert_i = points[i];
					int ttt = pointsIndex[i];
					points_New.push_back(points[i]);

					if (point_insert_i.size()!=3) {
						cout << "Hello";					
					}

					if (pointsIndex[i] != points_New.size() - 1) {

						cout << "Error index!";	

					}	

					vector<int> pointNeighbor_New_i;

					for (int j = 0; j < pointNeighbor[i].size()/2; j++) {

						int b1 = pointNeighbor[i][2 * j];
						int b2 = pointNeighbor[i][2 * j + 1];
						if (b1 == -1 && b2 == -1) {
							continue;						
						}
						else if (pointsIndex[b1] == -1 && pointsIndex[b2] == -1) {
							continue;
						}
						else if (b1 == -1 || b2 == -1 || pointsIndex[b1] == -1 || pointsIndex[b2] == -1) {
							cout << "Error, the index label is wrong!";						
						}
						else {							
							pointNeighbor_New_i.push_back(pointsIndex[b1]);
							pointNeighbor_New_i.push_back(pointsIndex[b2]);						
						}					
					
					}
					pointNeighbor_New.push_back(pointNeighbor_New_i);
				}
			}
			points.clear();
			points = points_New;
			pointNeighbor.clear();
			pointNeighbor = pointNeighbor_New;			

			//5. Check point neighbor structure
			for (int i = 0; i < pointNeighbor.size(); i++) {
				int b1 = i;
				vector<int> pointNeighbor_i = pointNeighbor[i];
				for (int j = 0; j < pointNeighbor_i.size() / 2; j++) {

					int b2 = pointNeighbor_i[2 * j];
					int b3 = pointNeighbor_i[2 * j + 1];

					if (b2 > pointNeighbor.size() || b2 < 0) {
						cout << "Hello";					
					}
					if (b3 > pointNeighbor.size() || b3 < 0) {
						cout << "Hello";
					}

					vector<int> pointNeighbor_i2 = pointNeighbor[b2];
					vector<int> pointNeighbor_i3 = pointNeighbor[b3];
					bool judge_b2 = false;
					bool judge_b3 = false;

					for (int k = 0; k < pointNeighbor_i2.size() / 2; k++) {

						int b21 = pointNeighbor_i2[2 * k];
						int b22 = pointNeighbor_i2[2 * k + 1];
						if (b21 == b1 && b22 == b3 || b21 == b3 && b22 == b1) {
							judge_b2 = true;
							break;
						}
					
					}

					for (int k = 0; k < pointNeighbor_i3.size() / 2; k++) {

						int b31 = pointNeighbor_i3[2 * k];
						int b32 = pointNeighbor_i3[2 * k + 1];
						if (b31 == b1 && b32 == b2 || b31 == b2 && b32 == b1) {
							judge_b3 = true;
							break;
						}

					}

					if (!(judge_b2 || judge_b2)) {
						cout << "wrong trangular:" << b1 << "," << b2 << "," << b3 << endl;					
					}				
				}			
			}

			for (int i = 0; i < points.size(); i++) {
				if (points[i].size() != 3) {
					cout << "Hello";				
				}			
			}

			if (countWhile == 0) {
				break;
			}
			else {
				countWhile--;
			}

		}
	}

	void MeshOptimization_OM_Flip() {
		
		while (1) {

			//1. Check collapse border
			int flipNum = 0;
			vector<bool> pointJudgeActivate(points.size(), true);//record point can be processed in parallel structure
			vector<int> pointFlipParallel;
			for (int i = 0; i < points.size(); i++) {
				if (pointJudgeActivate[i]) {
					pointFlipParallel.push_back(i);
					vector<int> points_index = pointNeighbor[i];
					for (int j = 0; j < points_index.size(); j++) {
						pointJudgeActivate[points_index[j]] = false;
					}
				}
			}

#pragma omp parallel for
			for (int i = 0; i < pointFlipParallel.size(); i++) {

				//Flip test					
				int b1 = pointFlipParallel[i];
				//cout << b1;			
				//Flip start
				vector<int> pNb1;//add points into list
				for (int j = 0; j < pointNeighbor[b1].size(); j++) {
					int b2_j = pointNeighbor[b1][j];
					if (!MeshOptimization_Exist(b2_j, pNb1)) {
						pNb1.push_back(b2_j);
					}
				}

				for (int j = 0; j < pNb1.size(); j++) {
					int b2 = pNb1[j];
					if (b2 < b1) {
						continue;
					}
					//cout << b2;
					vector<int> b3List;//store the other two points for flip
					for (int k = 0; k < pointNeighbor[b1].size() / 2; k++) {
						int b11 = pointNeighbor[b1][2 * k];
						int b12 = pointNeighbor[b1][2 * k + 1];
						if (b11 == b2 && !MeshOptimization_Exist(b12, b3List)) {
							b3List.push_back(b12);
						}
						if (b12 == b2 && !MeshOptimization_Exist(b11, b3List)) {
							b3List.push_back(b11);
						}
					}
					if (b3List.size() == 2) {
						//Flip processing
						//compute 4 point valence
						int b3 = b3List[0];
						int b4 = b3List[1];
						//cout <<"Flip:"<< b1 << "," << b2 << "," << b3 << "," << b4 << endl;
						if (MeshOptimization_If_Flip(b1, b2, b3, b4)) {
							//Flip update
							//cout << b1 << "," << b2 << "," << b3 << "," << b4 << endl;
							MeshOptimization_Flip_Processing(b1, b2, b3, b4);
							flipNum++;
						}
					}
					else {
						continue;
					}
				}			
			}
		    
			cout << flipNum;
			if (flipNum == 0) {
				cout << ".";
				break;
			}
			else {
				cout << ",";
			}		
		}
		cout << endl;

	}

	void MeshOptimization_OM_TangentSmoothing() {

		double landa = 0.6;
		vector<vector<double>> pointNew = points;

#pragma omp parallel for
		for (int i = 0; i < points.size(); i++) {
			vector<double> pi = points[i];
			int p1 = i;
			if (MeshOptimization_Exist(p1, points_Keep)) {
				//pointNew[i] = pi;
				continue;
			}
			vector<int> pointsN = pointNeighbor[i];
			vector<int> pointsN_withoutzRepeat;//neighbor points without repeat point
			vector<double> pointN_Area;//correspoinding area of neighbor point
			vector<double> pointB_List;//Adaptive L ditance for neighbor
			double AreaSum = 0;
			for (int j = 0; j < pointsN.size() / 2; j++) {
				int p2 = pointsN[2 * j];
				int p3 = pointsN[2 * j + 1];
				double areai = MeshOptimization_TrangularArea(p1, p2, p3);
				//double areai2 = areai / (double)2;
				int indexp2 = MeshOptimization_Exist_Index(p2, pointsN_withoutzRepeat);
				int indexp3 = MeshOptimization_Exist_Index(p3, pointsN_withoutzRepeat);
				if (indexp2 == -1) {
					pointsN_withoutzRepeat.push_back(p2);
					pointN_Area.push_back(areai);
				}
				else {
					pointN_Area[indexp2] = (pointN_Area[indexp2] + areai) / 2;
				}
				if (indexp3 == -1) {
					pointsN_withoutzRepeat.push_back(p3);
					pointN_Area.push_back(areai);
				}
				else {
					pointN_Area[indexp3] = (pointN_Area[indexp3] + areai) / 2;
				}
			}

			for (int j = 0; j < pointsN_withoutzRepeat.size(); j++) {

				double L_bi = MeshOptimization_KdTreeSearch(pointsN_withoutzRepeat[j]);
				pointB_List.push_back(L_bi);

			}

			for (int j = 0; j < pointN_Area.size(); j++) {

				AreaSum = AreaSum + pointN_Area[j] * pointB_List[j];

			}

			vector<double> gi(3, 0);
			for (int j = 0; j < pointsN_withoutzRepeat.size(); j++) {
				gi[0] = gi[0] + ((pointN_Area[j] * pointB_List[j]) / AreaSum) * points[pointsN_withoutzRepeat[j]][0];
				gi[1] = gi[1] + ((pointN_Area[j] * pointB_List[j]) / AreaSum) * points[pointsN_withoutzRepeat[j]][1];
				gi[2] = gi[2] + ((pointN_Area[j] * pointB_List[j]) / AreaSum) * points[pointsN_withoutzRepeat[j]][2];
			}

			//achieve normal of p1
			vector<double> ni = MeshOptimization_Normal_Point(p1);
			//achieve move vertor v_pg for gi and pi
			vector<double> v_pg(3);
			v_pg[0] = gi[0] - pi[0];
			v_pg[1] = gi[1] - pi[1];
			v_pg[2] = gi[2] - pi[2];
			double mapN_dis = v_pg[0] * ni[0] + v_pg[1] * ni[1] + v_pg[2] * ni[2];//v_pg map into n distance.
			vector<double> pi_new(3, 0);
			pi_new[0] = pi[0] + landa * (v_pg[0] - mapN_dis * ni[0]);
			pi_new[1] = pi[1] + landa * (v_pg[1] - mapN_dis * ni[1]);
			pi_new[2] = pi[2] + landa * (v_pg[2] - mapN_dis * ni[2]);
			if (!(pi_new[0]<9999 && pi_new[0]>-9999)) {
				continue;
			}
			pointNew[i] = pi_new;
		}

#pragma omp parallel for
		for (int i = 0; i < points.size(); i++) {
			points[i][0] = pointNew[i][0];
			points[i][1] = pointNew[i][1];
			points[i][2] = pointNew[i][2];
		}

	}

	//achieve the trangualr area
	double MeshOptimization_TrangularArea(int b1, int b2, int b3) {
		vector<double> p1 = points[b1];
		vector<double> p2 = points[b2];
		vector<double> p3 = points[b3];

		vector<double> v1(3);
		v1[0] = p2[0] - p1[0];
		v1[1] = p2[1] - p1[1];
		v1[2] = p2[2] - p1[2];

		vector<double> v2(3);
		v2[0] = p3[0] - p2[0];
		v2[1] = p3[1] - p2[1];
		v2[2] = p3[2] - p2[2];

		vector<double> v3(3);
		v3[0] = p1[0] - p3[0];
		v3[1] = p1[1] - p3[1];
		v3[2] = p1[2] - p3[2];

		double v1d = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
		double v2d = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
		double v3d = sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);

		double p = (v1d + v2d + v3d) / 2;
		double S = sqrt(p * (p - v1d) * (p - v2d) * (p - v3d));

		return S;
	}

	vector<double> MeshOptimization_Normal_Point(int b1) {

		vector<int> p1N = pointNeighbor[b1];
		vector<double> p1N_Weight;//weights of the normal list
		vector<vector<double>> p1N_Normal;//normal list
		double Area_sum = 0;
		for (int i = 0; i < p1N.size() / 2; i++) {
			int b2 = p1N[2 * i];
			int b3 = p1N[2 * i + 1];
			vector<double> n_i = MeshOptimization_Normal_Face(b1, b2, b3);
			double area_i = MeshOptimization_TrangularArea(b1, b2, b3);
			Area_sum = Area_sum + area_i;
			p1N_Weight.push_back(area_i);
			p1N_Normal.push_back(n_i);
		}
		vector<double> n(3, 0);
		for (int i = 0; i < p1N_Weight.size(); i++) {
			n[0] = n[0] + (p1N_Weight[i] / Area_sum) * p1N_Normal[i][0];
			n[1] = n[1] + (p1N_Weight[i] / Area_sum) * p1N_Normal[i][1];
			n[2] = n[2] + (p1N_Weight[i] / Area_sum) * p1N_Normal[i][2];
		}
		return n;
	}

	vector<double> MeshOptimization_Normal_Face(int b1, int b2, int b3) {

		vector<double> p1 = points[b1];
		vector<double> p2 = points[b2];
		vector<double> p3 = points[b3];

		vector<double> v1(3, 0);
		v1[0] = p2[0] - p1[0];
		v1[1] = p2[1] - p1[1];
		v1[2] = p2[2] - p1[2];

		vector<double> v2(3, 0);
		v2[0] = p3[0] - p2[0];
		v2[1] = p3[1] - p2[1];
		v2[2] = p3[2] - p2[2];

		vector<double> n(3, 0);
		n[0] = v2[2] * v1[1] - v2[1] * v1[2];
		n[1] = -v2[2] * v1[0] + v2[0] * v1[2];
		n[2] = v2[1] * v1[0] - v2[0] * v1[1];

		vector<double> n_Unit = MeshOptimization_Unit_Normal(n);
		return n_Unit;

	}

	vector<double> MeshOptimization_Unit_Normal(vector<double> n) {

		double length = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		n[0] = n[0] / length;
		n[1] = n[1] / length;
		n[2] = n[2] / length;
		return n;

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

	//Flip judge and processing
	bool MeshOptimization_If_Flip(int b1, int b2, int b3, int b4) {

		double disb34 = MeshOptimization_EdgeLength(b3, b4);

		double L_j1 = MeshOptimization_KdTreeSearch(points[b3]);
		double L_j2 = MeshOptimization_KdTreeSearch(points[b4]);
		if (L_j1 > L_j2) {
			L_j1 = L_j2;
		}
		if (disb34 > 1.33 * L_j1 || disb34 < 0.8 * L_j1) {
			return false;
		}
		bool b1e = MeshOptimization_Edge_Point(b1);
		bool b2e = MeshOptimization_Edge_Point(b2);
		bool b3e = MeshOptimization_Edge_Point(b3);
		bool b4e = MeshOptimization_Edge_Point(b4);
		int b1ValenceTarget = 6;
		if (b1e) {
			b1ValenceTarget = 4;
		}
		int b2ValenceTarget = 6;
		if (b2e) {
			b2ValenceTarget = 4;
		}
		int b3ValenceTarget = 6;
		if (b3e) {
			b3ValenceTarget = 4;
		}
		int b4ValenceTarget = 6;
		if (b4e) {
			b4ValenceTarget = 4;
		}
		int b1ValenceBefore = MeshOptimization_Point_Valence(b1);
		int b2ValenceBefore = MeshOptimization_Point_Valence(b2);
		int b3ValenceBefore = MeshOptimization_Point_Valence(b3);
		int b4ValenceBefore = MeshOptimization_Point_Valence(b4);
		int b1ValenceAfter = b1ValenceBefore - 1;
		int b2ValenceAfter = b2ValenceBefore - 1;
		int b3ValenceAfter = b3ValenceBefore + 1;
		int b4ValenceAfter = b4ValenceBefore + 1;

		int BeforeValence = 0;
		BeforeValence = abs(b1ValenceBefore - b1ValenceTarget) +
			abs(b2ValenceBefore - b2ValenceTarget) +
			abs(b3ValenceBefore - b3ValenceTarget) +
			abs(b4ValenceBefore - b4ValenceTarget);

		int AfterValence = 0;
		AfterValence = abs(b1ValenceAfter - b1ValenceTarget) +
			abs(b2ValenceAfter - b2ValenceTarget) +
			abs(b3ValenceAfter - b3ValenceTarget) +
			abs(b4ValenceAfter - b4ValenceTarget);
		//object Valence	

		if (AfterValence < BeforeValence) {
			return true;
		}
		else {
			return false;
		}
	}

	void MeshOptimization_Flip_Processing(int b1, int b2, int b3, int b4) {

		vector<int> b1N = pointNeighbor[b1];
		vector<int> b2N = pointNeighbor[b2];
		vector<int> b3N = pointNeighbor[b3];
		vector<int> b4N = pointNeighbor[b4];

		//update b1 structure
		vector<int> b1N_New;
		for (int i = 0; i < b1N.size() / 2; i++) {
			int b1N2 = b1N[2 * i];
			int b1N3 = b1N[2 * i + 1];
			if (b1N2 == b2 && (b1N3 == b3 || b1N3 == b4)) {
				b1N[2 * i] = -1;
				b1N[2 * i + 1] = -1;
				if (b1N3 == b3) {
					b1N_New.push_back(b4);
					b1N_New.push_back(b3);
				}
				else {
					b1N_New.push_back(b3);
					b1N_New.push_back(b4);
				}
			}
			if (b1N3 == b2 && (b1N2 == b3 || b1N2 == b4)) {
				b1N[2 * i] = -1;
				b1N[2 * i + 1] = -1;
				if (b1N2 == b3) {
					b1N_New.push_back(b3);
					b1N_New.push_back(b4);
				}
				else {
					b1N_New.push_back(b4);
					b1N_New.push_back(b3);
				}
			}
		}
		if (b1N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b1N_final;
		for (int i = 0; i < b1N.size() / 2; i++) {
			int b1N2 = b1N[2 * i];
			int b1N3 = b1N[2 * i + 1];
			if (b1N2 != -1 && b1N3 != -1) {
				b1N_final.push_back(b1N2);
				b1N_final.push_back(b1N3);
			}
		}
		b1N_final.push_back(b1N_New[0]);
		b1N_final.push_back(b1N_New[1]);

		//update b2 structure
		vector<int> b2N_New;
		for (int i = 0; i < b2N.size() / 2; i++) {
			int b2N2 = b2N[2 * i];
			int b2N3 = b2N[2 * i + 1];
			if (b2N2 == b1 && (b2N3 == b3 || b2N3 == b4)) {
				b2N[2 * i] = -1;
				b2N[2 * i + 1] = -1;
				if (b2N3 == b3) {
					b2N_New.push_back(b4);
					b2N_New.push_back(b3);
				}
				else {
					b2N_New.push_back(b3);
					b2N_New.push_back(b4);
				}
			}
			if (b2N3 == b1 && (b2N2 == b3 || b2N2 == b4)) {
				b2N[2 * i] = -1;
				b2N[2 * i + 1] = -1;
				if (b2N2 == b3) {
					b2N_New.push_back(b3);
					b2N_New.push_back(b4);
				}
				else {
					b2N_New.push_back(b4);
					b2N_New.push_back(b3);
				}
			}
		}
		if (b2N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b2N_final;
		for (int i = 0; i < b2N.size() / 2; i++) {
			int b2N2 = b2N[2 * i];
			int b2N3 = b2N[2 * i + 1];
			if (b2N2 != -1 && b2N3 != -1) {
				b2N_final.push_back(b2N2);
				b2N_final.push_back(b2N3);
			}
		}
		b2N_final.push_back(b2N_New[0]);
		b2N_final.push_back(b2N_New[1]);

		//update b3 structure	
		vector<int> b3N_New;
		for (int i = 0; i < b3N.size() / 2; i++) {
			int b3N2 = b3N[2 * i];
			int b3N3 = b3N[2 * i + 1];
			if ((b3N2 == b1 || b3N2 == b2) && (b3N3 == b1 || b3N3 == b2)) {
				b3N[2 * i] = -1;
				b3N[2 * i + 1] = -1;
				b3N_New.push_back(b3N2);
				b3N_New.push_back(b4);
				b3N_New.push_back(b4);
				b3N_New.push_back(b3N3);
			}
		}
		if (b3N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b3N_final;
		for (int i = 0; i < b3N.size() / 2; i++) {
			int b3N2 = b3N[2 * i];
			int b3N3 = b3N[2 * i + 1];
			if (b3N2 != -1 && b3N3 != -1) {
				b3N_final.push_back(b3N2);
				b3N_final.push_back(b3N3);
			}
		}
		b3N_final.push_back(b3N_New[0]);
		b3N_final.push_back(b3N_New[1]);
		b3N_final.push_back(b3N_New[2]);
		b3N_final.push_back(b3N_New[3]);

		//update b4 structure	
		vector<int> b4N_New;
		for (int i = 0; i < b4N.size() / 2; i++) {
			int b4N2 = b4N[2 * i];
			int b4N3 = b4N[2 * i + 1];
			if ((b4N2 == b1 || b4N2 == b2) && (b4N3 == b1 || b4N3 == b2)) {
				b4N[2 * i] = -1;
				b4N[2 * i + 1] = -1;
				b4N_New.push_back(b4N2);
				b4N_New.push_back(b3);
				b4N_New.push_back(b3);
				b4N_New.push_back(b4N3);
			}
		}
		if (b4N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b4N_final;
		for (int i = 0; i < b4N.size() / 2; i++) {
			int b4N2 = b4N[2 * i];
			int b4N3 = b4N[2 * i + 1];
			if (b4N2 != -1 && b4N3 != -1) {
				b4N_final.push_back(b4N2);
				b4N_final.push_back(b4N3);
			}
		}
		b4N_final.push_back(b4N_New[0]);
		b4N_final.push_back(b4N_New[1]);
		b4N_final.push_back(b4N_New[2]);
		b4N_final.push_back(b4N_New[3]);

		pointNeighbor[b1].clear();
		pointNeighbor[b2].clear();
		pointNeighbor[b3].clear();
		pointNeighbor[b4].clear();

		pointNeighbor[b1] = b1N_final;
		pointNeighbor[b2] = b2N_final;
		pointNeighbor[b3] = b3N_final;
		pointNeighbor[b4] = b4N_final;

	}

	int MeshOptimization_Point_Valence(int p1) {

		vector<int> b_num;
		for (int i = 0; i < pointNeighbor[p1].size() / 2; i++) {
			int p2 = pointNeighbor[p1][2 * i];
			int p3 = pointNeighbor[p1][2 * i + 1];
			if (!MeshOptimization_Exist(p2, b_num)) {
				b_num.push_back(p2);
			}
			if (!MeshOptimization_Exist(p3, b_num)) {
				b_num.push_back(p3);
			}
		}
		return b_num.size();

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
	
	void  AnIsotropic_OM_CheckErrorBorder() {//delete a border which is included in three faces at same time

		vector<int> pointActive(points.size());
		for (int i = 0; i < points.size(); i++) {
			pointActive[i] = i;		
		}
		for (int i = 0; i < pointNeighbor.size(); i++) {			
			vector<int> pointNeighbor_i = pointNeighbor[i];
			if (pointNeighbor_i.size() <= 2) {

				int b1 = pointNeighbor_i[0];
				int b2 = pointNeighbor_i[1];
				vector<int> b1n = pointNeighbor[b1];
				int count = 0;
				for (int j = 0; j < b1n.size(); j++) {
					if (b1n[j] == b2) {
						count++;					
					}
				}
				if (count >= 3) {

					pointActive[i] = -1;
					for (int j = i + 1; j < points.size(); j++) {

						pointActive[j]--;

					}
					for (int j = 0; j < pointNeighbor[b1].size()/2; j++) {

						int b11 = pointNeighbor[b1][2 * j];
						int b12 = pointNeighbor[b1][2 * j + 1];
						if (b11 == i || b12 == i) {
							pointNeighbor[b1][2 * j] = -1;
							pointNeighbor[b1][2 * j + 1] = -1;
						
						}
					
					}
					for (int j = 0; j < pointNeighbor[b2].size()/2; j++) {

						int b21 = pointNeighbor[b2][2 * j];
						int b22 = pointNeighbor[b2][2 * j + 1];

						if (b21 == i || b22 == i) {
							pointNeighbor[b2][2 * j] = -1;
							pointNeighbor[b2][2 * j + 1] = -1;

						}

					}
				
				}
				
			}
			else if (pointNeighbor_i.size() == 4) {
				int b1 = pointNeighbor_i[0];
				int b2 = pointNeighbor_i[1];
				int b3 = pointNeighbor_i[2];
				int b4 = pointNeighbor_i[3];
				if ((b1 == b3 || b1 == b4) && (b2 == b3 || b2 == b4)) {

					int b11 = pointNeighbor_i[0];
					int b12 = pointNeighbor_i[1];
					vector<int> b1n = pointNeighbor[b11];
					int count = 0;
					for (int j = 0; j < b1n.size(); j++) {
						if (b1n[j] == b12) {
							count++;
						}
					}
					if (count >= 3) {
						pointActive[i] = -1;
						for (int j = i + 1; j < points.size(); j++) {
							pointActive[j]--;
						}	
						for (int j = 0; j < pointNeighbor[b11].size() / 2; j++) {

							int b111 = pointNeighbor[b11][2 * j];
							int b112 = pointNeighbor[b11][2 * j + 1];
							if (b111 == i || b112 == i) {
								pointNeighbor[b11][2 * j] = -1;
								pointNeighbor[b11][2 * j + 1] = -1;

							}

						}
						for (int j = 0; j < pointNeighbor[b12].size() / 2; j++) {

							int b221 = pointNeighbor[b12][2 * j];
							int b222 = pointNeighbor[b12][2 * j + 1];

							if (b221 == i || b222 == i) {
								pointNeighbor[b12][2 * j] = -1;
								pointNeighbor[b12][2 * j + 1] = -1;

							}
						}
					}
				}
			}
			else {	
				continue;								
			}			
		}

		vector<vector<double>> points_new;
		vector<vector<int>> pointNeighbor_new;
		for (int i = 0; i < points.size(); i++) {
			if (pointActive[i] < 0) {
				continue;			
			}
			else {
				points_new.push_back(points[i]);			
			}		
		}

		points.clear();
		points = points_new;

		for (int i = 0; i < pointNeighbor.size(); i++) {
			//cout << i<< ",";
			if (pointActive[i] < 0) {
				continue;			
			}
			else {
				vector<int> pointNeighbor_i;
				for (int j = 0; j < pointNeighbor[i].size()/2; j++) {
					int b1 = pointNeighbor[i][2 * j];
					int b2 = pointNeighbor[i][2 * j + 1];
					if (b1 < 0 || b2 < 0) {
						continue;					
					}
					else if (pointActive[b1] < 0 || pointActive[b2] < 0) {
						continue;					
					}
					else {
						pointNeighbor_i.push_back(pointActive[b1]);
						pointNeighbor_i.push_back(pointActive[b2]);
					}					
				}
				pointNeighbor_new.push_back(pointNeighbor_i);			
			}		
		}	

		pointNeighbor.clear();
		pointNeighbor = pointNeighbor_new;

		cout << endl;

	}

};


