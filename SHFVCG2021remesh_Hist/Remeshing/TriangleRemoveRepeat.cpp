/*********************************************************************************

							 Triangluar Point Check

							 Updating in 2021/05/21

							   By Dr. Chenlei Lv

			The functions includes:
			1. Remove repeat point in mesh
			2. Update mesh

**********************************************************************************/
#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm> 
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
using namespace std;

class Triangluar_RemoveRepeat {

private:

	vector<vector<double>> pointSet;
	vector<vector<int>> faceSet;
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;//kdtree
	int K = 10;

public:

	void TR_Init(vector<vector<double>> pointSetInput, vector<vector<int>> faceSetInput) {

		std::cout << "Triangluar_RemoveRepeat Init" << std::endl;
		//init data:

		pointSet = pointSetInput;
		faceSet = faceSetInput;

		//init kdtree:
		std::cout << "Init kdtree" << std::endl;
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->width = pointSet.size();
		cloud->height = 1;
		cloud->points.resize(cloud->width * cloud->height);
		// fills a PointCloud with random data
		for (int i = 0; i < pointSet.size(); i++)
		{
			pcl::PointXYZ pxyz;
			cloud->points[i].x = pointSet[i][0];
			cloud->points[i].y = pointSet[i][1];
			cloud->points[i].z = pointSet[i][2];

		}
		kdtree.setInputCloud(cloud);

	}

	void TR_Start(string filename) {

		std::cout << "Remove repeat points" << std::endl;
		TR_Detect();
		TR_RemoveFace();
		string filenamePath = "Data//subRemesh//" + filename + "_removeRepeat.obj";

		//TR_SaveOBJ(filenamePath, pointSet, faceSet);
		std::cout << "Triangluar_RemoveRepeat finished!" << std::endl;

	}

	vector<vector<double>> TR_Get_PointSet() {

		return pointSet;

	}

	vector<vector<int>> TR_Get_FaceSet() {

		return faceSet;

	}

private:

	void TR_Detect() {

		vector<int> pointIndex(pointSet.size());
		for (int i = 0; i < pointIndex.size(); i++) {
			pointIndex[i] = i;
		}

		for (int i = 0; i < pointSet.size(); i++) {
			int b_index = i;
			//if (b_index == 15689) {
				//cout << "Hello!" << endl;			
			//}
			if (pointIndex[b_index] != b_index) {
				continue;
			}
			std::vector<int> pointIdxNKNSearch(K);
			std::vector<float> pointNKNSquaredDistance(K);
			std::vector<int> pointNeibor_i;
			pcl::PointXYZ searchPoint;
			searchPoint.x = pointSet[i][0];
			searchPoint.y = pointSet[i][1];
			searchPoint.z = pointSet[i][2];
			kdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);

			vector<int> detectIndex;
			for (int j = 0; j < pointIdxNKNSearch.size(); j++) {
				if (pointNKNSquaredDistance[j] < 0.0001 && pointIdxNKNSearch[j] != b_index) {
					detectIndex.push_back(pointIdxNKNSearch[j]);
				}
			}

			if (detectIndex.size() >= 9) {
				detectIndex.clear();
				detectIndex = TR_Adaptive_K(b_index);
			}

			for (int j = 0; j < detectIndex.size(); j++) {
				pointIndex[detectIndex[j]] = b_index;
			}
		}

		//TR_Update(pointIndex);
		TR_Update2(pointIndex);
	}

	vector<int> TR_Adaptive_K(int b_Index) {

		int K_adaptive = 100;
		while (true) {
			std::vector<int> pointIdxNKNSearch(K_adaptive);
			std::vector<float> pointNKNSquaredDistance(K_adaptive);
			std::vector<int> pointNeibor_i;
			pcl::PointXYZ searchPoint;
			searchPoint.x = pointSet[b_Index][0];
			searchPoint.y = pointSet[b_Index][1];
			searchPoint.z = pointSet[b_Index][2];
			kdtree.nearestKSearch(searchPoint, K_adaptive, pointIdxNKNSearch, pointNKNSquaredDistance);
			vector<int> detectIndex;
			for (int j = 0; j < pointIdxNKNSearch.size(); j++) {
				if (pointNKNSquaredDistance[j] < 0.0001 && pointIdxNKNSearch[j] != b_Index) {
					detectIndex.push_back(pointIdxNKNSearch[j]);
				}
			}
			if (detectIndex.size() >= K_adaptive - 1) {
				K_adaptive = K_adaptive + 100;
			}
			else {
				return detectIndex;
			}
		}
	}


	void TR_Update(vector<int> pointIndex) {

		vector<vector<double>> pointSet_Update;


		int count = 0;
		for (int i = 0; i < pointIndex.size(); i++) {
			if (pointIndex[i] != i) {
				pointIndex[i] = pointIndex[pointIndex[i]];
				count++;
			}
			else {
				pointIndex[i] = pointIndex[i] - count;
				pointSet_Update.push_back(pointSet[i]);
			}
		}

		pointSet.clear();
		pointSet = pointSet_Update;

		for (int i = 0; i < faceSet.size(); i++) {

			int b1 = faceSet[i][0];
			int b2 = faceSet[i][1];
			int b3 = faceSet[i][2];
			faceSet[i][0] = pointIndex[b1];
			faceSet[i][1] = pointIndex[b2];
			faceSet[i][2] = pointIndex[b3];

		}

	}

	void TR_Update2(vector<int> pointIndex) {

		vector<vector<double>> pointSet_Update;
		vector<bool> pointIndexJudge(pointIndex.size(), true);

		for (int i = 0; i < pointIndex.size(); i++) {
			if (pointIndex[i] != i) {
				pointIndexJudge[i] = false;
			}
		}

		int count = 0;
		for (int i = 0; i < pointIndexJudge.size(); i++) {
			if (pointIndexJudge[i]) {
				pointSet_Update.push_back(pointSet[i]);
				pointIndex[i] = count;
				count++;
			}
		}

		for (int i = 0; i < pointIndexJudge.size(); i++) {
			if (!pointIndexJudge[i]) {
				int pIndex_i = pointIndex[i];
				pointIndex[i] = pointIndex[pIndex_i];
			}
		}

		pointSet.clear();
		pointSet = pointSet_Update;

		for (int i = 0; i < faceSet.size(); i++) {

			int b1 = faceSet[i][0];
			int b2 = faceSet[i][1];
			int b3 = faceSet[i][2];
			faceSet[i][0] = pointIndex[b1];
			faceSet[i][1] = pointIndex[b2];
			faceSet[i][2] = pointIndex[b3];

		}

	}

	//remove surface: three points in a line
	void TR_RemoveFace() {

		vector<vector<int>> faceSetNew;
		for (int i = 0; i < faceSet.size(); i++) {

			//if (i == 25131) {
				//cout << "i:" << i << endl;			
			//}
			int b1 = faceSet[i][0];
			int b2 = faceSet[i][1];
			int b3 = faceSet[i][2];

			//if (b1 == 58345 && b2 == 32264 && b3 == 45499) {
				//cout << "Hello!" << endl;			
			//}

			vector<double> v1(3);
			v1[0] = pointSet[b1][0] - pointSet[b2][0];
			v1[1] = pointSet[b1][1] - pointSet[b2][1];
			v1[2] = pointSet[b1][2] - pointSet[b2][2];

			vector<double> v2(3);
			v2[0] = pointSet[b2][0] - pointSet[b3][0];
			v2[1] = pointSet[b2][1] - pointSet[b3][1];
			v2[2] = pointSet[b2][2] - pointSet[b3][2];

			double angle_i = TC_VectorAngel(v1, v2);

			if (angle_i < 0.0005 || angle_i > 3.141) {
				continue;
			}
			else if (b1 == b2 || b1 == b3 || b2 == b3) {
				continue;
			}
			else {
				faceSetNew.push_back(faceSet[i]);
			}
		}
		faceSet.clear();
		faceSet = faceSetNew;
	}

	double TC_VectorAngel(vector<double> v1, vector<double> v2) {

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
		return acos(cosV12);

	}

	void TR_SaveOBJ(string fileName, vector<vector<double>> points, vector<vector<int>> facet) {

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
