/*********************************************************************************

				  Mesh Geometric Features Computation

						Updating in 2020/12/01

						   By Dr. Chenlei Lv

			The functions includes:
			1. Guassian Curvature
			2. Mean Curvature
			3. Laplace-B operator

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
#define PI 3.1415926
using namespace std;

class MeshGeometric {

private:

	int step = 3;
	double landa = 0.5;
	vector<vector<double>> points;
	vector<vector<int>> pNeighbor;
	vector<double> G_Value;
	vector<double> M_Value;
	vector<double> N_Value;
	vector<double> HN_Value;

public:

	void MeshGeometric_init(vector<vector<double>> points_input, vector<vector<int>> pNeighbor_input) {

		points = points_input;
		pNeighbor = pNeighbor_input;
		G_Value.resize(points.size());
		M_Value.resize(points.size());
		N_Value.resize(points.size());
		MeshGeometric_G_M_N_Value();
	}

	void MeshGeometric_init_Harmonic(vector<vector<double>> points_input, vector<vector<int>> pNeighbor_input) {

		points = points_input;
		pNeighbor = pNeighbor_input;
		N_Value.resize(points.size());
		HN_Value.resize(points.size());
		MeshGeometric_Harmonic_N_Value();
	}

	vector<double> MeshGeometric_Get_G() {

		return G_Value;

	}

	vector<double> MeshGeometric_Get_M() {

		return M_Value;

	}

	vector<double> MeshGeometric_Get_N() {

		return N_Value;

	}

	vector<double> MeshGeometric_Get_HN() {

		return HN_Value;

	}

private:

	void MeshGeometric_G_M_N_Value() {

		for (int i = 0; i < points.size(); i++) {
			//cout << "point:" << i << endl;
			//if (i == 1543) {
				//cout << "1543,hello!" << endl;			
			//}
			double innerAngleSum = 0;
			vector<int> p1_Neighbor = pNeighbor[i];
			vector<int> pointsN_withoutzRepeat;
			for (int j = 0; j < p1_Neighbor.size() / 2; j++) {
				int p2 = p1_Neighbor[2 * j];
				int p3 = p1_Neighbor[2 * j + 1];
				innerAngleSum = innerAngleSum + MeshGeometric_InnerAngle(i, p2, p3);
				//double areai2 = areai / (double)2;
				int indexp2 = MeshGeometric_Exist_Index(p2, pointsN_withoutzRepeat);
				int indexp3 = MeshGeometric_Exist_Index(p3, pointsN_withoutzRepeat);
				if (indexp2 == -1) {
					pointsN_withoutzRepeat.push_back(p2);
				}
				if (indexp3 == -1) {
					pointsN_withoutzRepeat.push_back(p3);
				}
			}
			double Gi = 2 * PI - innerAngleSum;
			G_Value[i] = Gi;

			vector<double> ctan_weight(pointsN_withoutzRepeat.size(), 0);
			for (int j = 0; j < pointsN_withoutzRepeat.size(); j++) {
				int p2 = pointsN_withoutzRepeat[j];
				//achieve related two points
				for (int k = 0; k < p1_Neighbor.size() / 2; k++) {

					int p21 = p1_Neighbor[2 * k];
					int p22 = p1_Neighbor[2 * k + 1];
					int p2_real;
					if (p21 != p2 && p22 != p2) {
						continue;
					}
					else {
						if (p21 == p2) {
							p2_real = p22;
						}
						else {
							p2_real = p21;
						}
						double angle = MeshGeometric_InnerAngle(p2_real, i, p2);
						double tanValue = tan(abs(angle));
						if (tanValue < 0.1) {
							tanValue = 0.1;
						}
						if (tanValue > 10) {
							tanValue = 10;
						}
						ctan_weight[j] = ctan_weight[j] + 1.0 / tan(angle);
					}
				}
			}
			//achieve cot curvature
			vector<double> vLB(3, 0);
			for (int j = 0; j < pointsN_withoutzRepeat.size(); j++) {
				vLB[0] = vLB[0] + ctan_weight[j] * (points[pointsN_withoutzRepeat[j]][0] - points[i][0]);
				vLB[1] = vLB[1] + ctan_weight[j] * (points[pointsN_withoutzRepeat[j]][1] - points[i][1]);
				vLB[2] = vLB[2] + ctan_weight[j] * (points[pointsN_withoutzRepeat[j]][2] - points[i][2]);
			}
			double Hi = sqrt(vLB[0] * vLB[0] + vLB[1] * vLB[1] + vLB[2] * vLB[2]) / 2;
			M_Value[i] = Hi;
		}

		for (int i = 0; i < points.size(); i++) {
			vector<double> points_i_n = MeshGeometric_Normal_Point(i);
			double angleSum = 0;
			for (int j = 0; j < pNeighbor[i].size() / 2; j++) {
				int p2 = pNeighbor[i][2 * j];
				int p3 = pNeighbor[i][2 * j + 1];
				vector<double> p2n = MeshGeometric_Normal_Point(p2);
				vector<double> p3n = MeshGeometric_Normal_Point(p3);
				double a2 = MeshGeometric_Angle(points_i_n, p2n);
				double a3 = MeshGeometric_Angle(points_i_n, p3n);
				angleSum = angleSum + a2 + a3;
			}
			N_Value[i] = angleSum / pNeighbor[i].size();
		}
	}

	void MeshGeometric_Harmonic_N_Value() {

		vector<vector<int>> pointsN_withoutzRepeat;
		vector<vector<double>> ctan_weight_v;
		vector<vector<double>> distance_v;

		//initial N and cotangent wij
		cout << "initial N and cotangent wij" << endl;
		for (int i = 0; i < points.size(); i++) {
			
			vector<double> points_i_n = MeshGeometric_Normal_Point(i);
			vector<double> points_i = points[i];
			vector<int> pointsN_withoutzRepeat_i;
			vector<double> distance_i;
			double angleSum = 0;
			for (int j = 0; j < pNeighbor[i].size() / 2; j++) {
				int p2 = pNeighbor[i][2 * j];
				int p3 = pNeighbor[i][2 * j + 1];
				vector<double> points_p2 = points[p2];
				vector<double> points_p3 = points[p3];
				double dpi2 = sqrt((points_i[0] - points_p2[0]) * (points_i[0] - points_p2[0])
					+ (points_i[1] - points_p2[1]) * (points_i[1] - points_p2[1])
					+ (points_i[2] - points_p2[2]) * (points_i[2] - points_p2[2]));
				double dpi3 = sqrt((points_i[0] - points_p3[0]) * (points_i[0] - points_p3[0])
					+ (points_i[1] - points_p3[1]) * (points_i[1] - points_p3[1])
					+ (points_i[2] - points_p3[2]) * (points_i[2] - points_p3[2]));

				vector<double> p2n = MeshGeometric_Normal_Point(p2);
				vector<double> p3n = MeshGeometric_Normal_Point(p3);

				if (p2n[0] == 0 && p2n[1] == 0 && p2n[2] == 0) {
					p2n.clear();
					p2n = p3n;				
				}

				if (p3n[0] == 0 && p3n[1] == 0 && p3n[2] == 0) {
					p3n.clear();
					p3n = p2n;				
				}				

				double a2 = MeshGeometric_Angle(points_i_n, p2n);
				double a3 = MeshGeometric_Angle(points_i_n, p3n);
				angleSum = angleSum + a2 + a3;
				int indexp2 = MeshGeometric_Exist_Index(p2, pointsN_withoutzRepeat_i);
				int indexp3 = MeshGeometric_Exist_Index(p3, pointsN_withoutzRepeat_i);
				if (indexp2 == -1) {
					pointsN_withoutzRepeat_i.push_back(p2);
					distance_i.push_back(dpi2);
				}
				if (indexp3 == -1) {
					pointsN_withoutzRepeat_i.push_back(p3);
					distance_i.push_back(dpi3);
				}
			}

			//if (!(angleSum < 10 && angleSum >-10)) {
				//cout << "Hello?";							
			//}

			vector<int> p1_Neighbor = pNeighbor[i];
			vector<double> ctan_weight(pointsN_withoutzRepeat_i.size(), 0);
			for (int j = 0; j < pointsN_withoutzRepeat_i.size(); j++) {
				//if (i == 967 && j == pointsN_withoutzRepeat_i.size()-1) {
					//cout << "Error!";
				//}
				int p2 = pointsN_withoutzRepeat_i[j];
				//achieve related two points
				for (int k = 0; k < p1_Neighbor.size() / 2; k++) {
					int p21 = p1_Neighbor[2 * k];
					int p22 = p1_Neighbor[2 * k + 1];
					int p2_real;
					if (p21 != p2 && p22 != p2) {
						continue;
					}
					else {
						if (p21 == p2) {
							p2_real = p22;
						}
						else {
							p2_real = p21;
						}
						double angle = MeshGeometric_InnerAngle(p2_real, i, p2);
						double tanValue = tan(abs(angle));
						if (tanValue < 0.1) {
							tanValue = 0.1;
						}
						if (tanValue > 10) {
							tanValue = 10;
						}
						ctan_weight[j] = ctan_weight[j] + 1.0 / tanValue;
					}
				}
			}
			ctan_weight_v.push_back(ctan_weight);
			pointsN_withoutzRepeat.push_back(pointsN_withoutzRepeat_i);
			distance_v.push_back(distance_i);
			N_Value[i] = angleSum / pNeighbor[i].size();
			HN_Value[i] = N_Value[i];

			//if (!(HN_Value[i] < 10 && HN_Value[i] >-10)) {
				//cout << "Hello?";
			//}
		}

		//initial cotangent wij with distanc weight
		cout << "initial cotangent wij with distanc weight" << endl;
		vector<vector<double>> ctan_weight_v_f;//ctan_weight with distance weighted.
		for (int i = 0; i < points.size(); i++) {			
			//cout << i << endl;
			vector<double> ctan_weight_i = ctan_weight_v[i];
			vector<double> ctan_weight_i_f(ctan_weight_i.size());
			vector<int> pointsN_withoutzRepeat_i = pointsN_withoutzRepeat[i];
			vector<double> distance_i = distance_v[i];
			double sum_weight = 0;
			for (int j = 0; j < pointsN_withoutzRepeat_i.size(); j++) {
				sum_weight = sum_weight + ctan_weight_i[j] / distance_i[j];
			}
			for (int j = 0; j < pointsN_withoutzRepeat_i.size(); j++) {
				ctan_weight_i_f[j] = (ctan_weight_i[j] / distance_i[j]) / sum_weight;
				//if (!(ctan_weight_i_f[j] < 10 && ctan_weight_i_f[j] >-10)) {
					//cout << "Hello?";
				//}
			}
			ctan_weight_v_f.push_back(ctan_weight_i_f);
			
		}

		//smooth processing
		cout << "smooth processing." << endl;
		while (step--) {
			for (int i = 0; i < points.size(); i++) {
				//if (i == 113) {
					//cout << "Error!";
				//}
				vector<double> ctan_weight_i = ctan_weight_v_f[i];
				vector<int> pointsN_withoutzRepeat_i = pointsN_withoutzRepeat[i];
				double HN_i = 0;
				for (int j = 0; j < pointsN_withoutzRepeat_i.size(); j++) {
					HN_i = HN_i + ctan_weight_i[j] * HN_Value[pointsN_withoutzRepeat_i[j]];
				}
				HN_Value[i] = HN_Value[i] * landa + HN_i * (1 - landa);
				//if (!(HN_Value[i] < 10 && HN_Value[i] >-10)) {
					//cout << "Hello?";
				//}
			}
		}

	}

	int MeshGeometric_Exist_Index(int b, vector<int> v) {

		for (int i = 0; i < v.size(); i++) {
			if (b == v[i]) {
				return i;
			}
		}
		return -1;
	}

	double MeshGeometric_InnerAngle(int p1, int p2, int p3) {

		vector<double> p1p = points[p1];
		vector<double> p2p = points[p2];
		vector<double> p3p = points[p3];

		vector<double> pn1(3);
		pn1[0] = p2p[0] - p1p[0];
		pn1[1] = p2p[1] - p1p[1];
		pn1[2] = p2p[2] - p1p[2];

		vector<double> pn2(3);
		pn2[0] = p3p[0] - p1p[0];
		pn2[1] = p3p[1] - p1p[1];
		pn2[2] = p3p[2] - p1p[2];

		double pn1l = sqrt(pn1[0] * pn1[0] + pn1[1] * pn1[1] + pn1[2] * pn1[2]);
		pn1[0] = pn1[0] / pn1l;
		pn1[1] = pn1[1] / pn1l;
		pn1[2] = pn1[2] / pn1l;

		double pn2l = sqrt(pn2[0] * pn2[0] + pn2[1] * pn2[1] + pn2[2] * pn2[2]);
		pn2[0] = pn2[0] / pn2l;
		pn2[1] = pn2[1] / pn2l;
		pn2[2] = pn2[2] / pn2l;

		double cosa1 = pn1[0] * pn2[0] + pn1[1] * pn2[1] + pn1[2] * pn2[2];
		if (cosa1 > 1) {
			cosa1 = 1;		
		}
		if (cosa1 < -1) {
			cosa1 = -1;		
		}
		double anglea1 = acos(cosa1);
		//if (!(anglea1 <= 3.1415 && anglea1 >= -3.1415)) {
			//cout << "Hello!";		
		//}
		return anglea1;


	}

	vector<double> MeshGeometric_Normal_Point(int b1) {

		vector<int> p1N = pNeighbor[b1];
		vector<double> p1N_Weight;//weights of the normal list
		vector<vector<double>> p1N_Normal;//normal list
		double Area_sum = 0;
		for (int i = 0; i < p1N.size() / 2; i++) {
			int b2 = p1N[2 * i];
			int b3 = p1N[2 * i + 1];
			vector<double> n_i = MeshGeometric_Normal_Face(b1, b2, b3);
			double area_i = MeshGeometric_TrangularArea(b1, b2, b3);
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

	vector<double> MeshGeometric_Normal_Face(int b1, int b2, int b3) {

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

		vector<double> n_Unit = MeshGeometric_Unit_Normal(n);
		return n_Unit;

	}

	vector<double> MeshGeometric_Unit_Normal(vector<double> n) {

		double length = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		if (length == 0) {
			return n;		
		}
		n[0] = n[0] / length;
		n[1] = n[1] / length;
		n[2] = n[2] / length;
		return n;

	}

	double MeshGeometric_TrangularArea(int b1, int b2, int b3) {
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

	double MeshGeometric_Angle(vector<double> v1, vector<double> v2) {

		double cosa1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
		if (cosa1 > 1) {
			cosa1 = 1;		
		}
		if (cosa1 < -1) {
			cosa1 = -1;		
		}
		double anglea1 = acos(cosa1);
		return anglea1;

	}

	bool MeshGeometric_Exist(int b, vector<int> v) {

		for (int i = 0; i < v.size(); i++) {
			if (b == v[i]) {
				return true;
			}
		}
		return false;
	}
};


