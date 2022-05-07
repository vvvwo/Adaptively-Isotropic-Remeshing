/*********************************************************************************

					 Geodesic Computation in a point cloud Pro
								by Fast Marching
			   (Improve the function for Non-continuous Point clouds)


							 Updating in 2021/06/24

							   By Dr. Chenlei Lv

			The functions includes:
			1. Using the Fast Marching to achieve the geodeisc map from
			one source (point index or insert a point)
			2. Using the Fast Marching to achieve the geodeisc map for
			Non-continuous Point clouds
			3. Rebuild the geodesic path between two points;


*********************************************************************************/

#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h> 
#include <time.h>
#include <functional>
#include <algorithm>
#include<vcg/complex/complex.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>

using namespace std;
using namespace vcg;

class pPairDic {//store neighbor structure
public:
	int source;
	int target;
	float distance;
};

class PointGeodesic_Pro {

private:

	vector<Point3f> pointSet;
	vector<Point3f> pointNormal;
	pcl::KdTreeFLANN<pcl::PointXYZ> pkdtree;
	vector<int> pointLabel;
	vector<float> pointDistance;
	int pIndexGlobal;
	vector<vector<pPairDic>> PointPairDictionary;


public:

	void PointGeodesic_Pro_init(vector<Point3f> pointSetInput, vector<Point3f> pointNormalInput) {

		pointSet = pointSetInput;//point position
		pointNormal = pointNormalInput;//point normal
		pkdtree = PointGeodesic_KdTree(pointSet);//contruct kdtree
		PointGeodesic_Init_PointPairDictionary(0);//init PointPairDictionary
	}

	vector<int> PointGeodesic_Pro_ReturnAccuratePoints(vector<Point3f> pointSetInput, vector<Point3f> pointNormalInput, int subNumber) {

		pointSet = pointSetInput;//point position
		pointNormal = pointNormalInput;//point normal
		pkdtree = PointGeodesic_KdTree(pointSet);//contruct kdtree

		//init global searching distance and pointLabel
		if (pointLabel.size() > 0 || pointDistance.size() > 0) {
			pointLabel.clear();
			pointDistance.clear();
		}

		if (PointPairDictionary.size() > 0) {
			PointPairDictionary.clear();
		}

		int pIndex = 0;
		pIndexGlobal = pIndex;
		pointLabel.resize(pointSet.size(), 2);//0:Alive; 1:Close; 2:Far 
		PointPairDictionary.resize(pointSet.size());
		pointDistance.resize(pointSet.size(), 99999);
		if (pIndex >= pointLabel.size() || pIndex < 0) {
			cout << "label input error, out of vector";
			vector<int> t;
			return t;
		}

		int iterationNum = pointSet.size();
		vector<vector<int>> subset;
		//Achieve differnet continuous point cloud set
		int stepNum = 1;
		while (true) {
			cout << "step:" << stepNum << endl;
			stepNum++;
			vector<int> subset_i = PointGeodesic_DisMap_Unit(pIndex, iterationNum);
			subset.push_back(subset_i);
			iterationNum = iterationNum - subset_i.size();
			if (iterationNum <= 0) {
				break;
			}
			else {
				bool subJudge = true;
				for (int i = 0; i < pointLabel.size(); i++) {
					if (pointLabel[i] == 2) {
						pIndex = i;
						subJudge = false;
					}
				}
				if (subJudge) {
					break;//all points are assigned value.				
				}
			}
			cout << endl;
		}
		int maxSub = 0;
		int maxSuNum = subset[0].size();
		for (int i = 0; i < subset.size(); i++) {
			if (subset[i].size() > maxSuNum) {
				maxSuNum = subset[i].size();
				maxSub = i;
			}
		}
		//assigh class index for all point
		vector<int> pointClass(pointSet.size(), -1);

		//update subset, delete small set of point cloud part;
		vector<vector<int>> subset_new;
		for (int i = 0; i < subset.size(); i++) {
			if (i == maxSub) {
				maxSub = subset_new.size();
			}
			if (subset[i].size() > subNumber) {
				int classIndex = subset_new.size();
				subset_new.push_back(subset[i]);
				for (int j = 0; j < subset[i].size(); j++) {
					int indexij = subset[i][j];
					if (indexij < 0 || indexij >= pointClass.size()) {
						cout << "Point Index Error" << endl;
					}
					else {
						pointClass[indexij] = classIndex;
					}
				}
			}
		}
		subset.clear();
		subset = subset_new;
		vector<int> resultFinal;
		for (int i = 0; i < subset.size(); i++) {
			for (int j = 0; j < subset[i].size(); j++) {
				int pIndex = subset[i][j];
				resultFinal.push_back(pIndex);
			}
		}
		return resultFinal;

	}

	vector<float> PointGeodesic_DisMap(int pIndex) {

		if (pointLabel.size() > 0 || pointDistance.size() > 0) {
			pointLabel.clear();
			pointDistance.clear();
		}

		pIndexGlobal = pIndex;
		pointLabel.resize(pointSet.size(), 2);//0:Alive; 1:Close; 2:Far 
		pointDistance.resize(pointSet.size(), 99999);
		if (pIndex >= pointLabel.size() || pIndex < 0) {
			cout << "label input error, out of vector";
			return pointDistance;
		}
		pointLabel[pIndex] = 0;
		pointDistance[pIndex] = 0;
		bool judge = true;

		//init:
		vector<pair<float, int>> pointCloseDistanceHeap;//Close Set min-heap
		int K = 7;
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);
		pcl::PointXYZ searchPoint;
		searchPoint.x = pointSet[pIndex][0];
		searchPoint.y = pointSet[pIndex][1];
		searchPoint.z = pointSet[pIndex][2];
		pkdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);

		for (int i = 0; i < PointPairDictionary[pIndex].size(); i++) {
			int pi_source = PointPairDictionary[pIndex][i].source;
			int pi_target = PointPairDictionary[pIndex][i].target;
			float pi_dis = PointPairDictionary[pIndex][i].distance;
			if (pi_source == pIndex) {
				pointIdxNKNSearch.push_back(pi_target);
				pointNKNSquaredDistance.push_back(pi_dis * pi_dis);
			}
			else if (pi_target == pIndex) {
				pointIdxNKNSearch.push_back(pi_source);
				pointNKNSquaredDistance.push_back(pi_dis * pi_dis);
			}
			else {
				cout << "PointPairDictionary add error!" << endl;
			}
		}


		for (int i = 1; i < pointIdxNKNSearch.size(); i++) {
			int index_i = pointIdxNKNSearch[i];
			pointLabel[index_i] = 1;
			double distance_i = sqrt(pointNKNSquaredDistance[i]);
			pointDistance[index_i] = distance_i;
			pointCloseDistanceHeap.push_back(make_pair(distance_i, index_i));
		}



		//init min-heap		
		make_heap(pointCloseDistanceHeap.begin(), pointCloseDistanceHeap.end(), greater<pair<float, int>>());

		int iterationNum = pointSet.size();
		int iterstep = iterationNum / 10;
		//Fast Mearching
		while (judge) {
			iterationNum--;
			//0:Alive; 1:Close; 2:Far 

			//achieve the min value of the min-heap
			pop_heap(pointCloseDistanceHeap.begin(), pointCloseDistanceHeap.end(), greater<pair<float, int>>());
			pair<float, int> t = pointCloseDistanceHeap[pointCloseDistanceHeap.size() - 1];
			int indexMin = t.second;
			//if (indexMin == 0) {
				//cout << "Hello!" << endl;
			//}
			pointLabel[indexMin] = 0;//Set the point into Alive
			pointCloseDistanceHeap.pop_back();

			//achive the new point's neighbor
			std::vector<int> pointIdxNKNSearch_t(K);
			std::vector<float> pointNKNSquaredDistance_t(K);
			pcl::PointXYZ searchPoint_t;
			searchPoint_t.x = pointSet[indexMin][0];
			searchPoint_t.y = pointSet[indexMin][1];
			searchPoint_t.z = pointSet[indexMin][2];
			pkdtree.nearestKSearch(searchPoint_t, K, pointIdxNKNSearch_t, pointNKNSquaredDistance_t);

			for (int i = 0; i < PointPairDictionary[indexMin].size(); i++) {
				int pi_source = PointPairDictionary[indexMin][i].source;
				int pi_target = PointPairDictionary[indexMin][i].target;
				float pi_dis = PointPairDictionary[indexMin][i].distance;
				if (pi_source == indexMin) {
					pointIdxNKNSearch_t.push_back(pi_target);
					pointNKNSquaredDistance_t.push_back(pi_dis * pi_dis);
				}
				else if (pi_target == indexMin) {
					pointIdxNKNSearch_t.push_back(pi_source);
					pointNKNSquaredDistance_t.push_back(pi_dis * pi_dis);
				}
				else {
					cout << "While: PointPairDictionary add error!" << endl;
				}
			}

			vector<int> activeSet;
			vector<float> activeSetDis;
			vector<int> closeSet;
			vector<float> closeSetDis;
			vector<int> farSet;
			vector<float> farSetDis;

			//update the neighbor distance 
			for (int i = 1; i < pointIdxNKNSearch_t.size(); i++) {
				int indexi_n = pointIdxNKNSearch_t[i];
				float dis_t = sqrt(pointNKNSquaredDistance_t[i]);

				if (pointLabel[indexi_n] == 1) {
					closeSet.push_back(indexi_n);
					float closeSetDis_i = PointGeodesic_DistanceSinglePoint(indexi_n);
					if (!(closeSetDis_i <= 99999 && closeSetDis_i >= 0)) {
						cout << "Hello!" << endl;
						//pointDistance[indexi_n] = 0;
					}

					if (closeSetDis_i < pointDistance[indexi_n]) {
						pointDistance[indexi_n] = closeSetDis_i;
					}
					else {
						closeSetDis_i = pointDistance[indexi_n];
					}
					closeSetDis.push_back(closeSetDis_i);
				}
				else if (pointLabel[indexi_n] == 0) {
					activeSet.push_back(indexi_n);
					activeSetDis.push_back(pointDistance[indexi_n]);
				}
				else {
					farSet.push_back(indexi_n);
					float farSetDis_i = PointGeodesic_DistanceSinglePoint(indexi_n);
					if (farSetDis_i > dis_t + pointDistance[indexMin]) {
						farSetDis_i = dis_t + pointDistance[indexMin];
					}
					farSetDis.push_back(farSetDis_i);
				}
			}

			float indexMinD = PointGeodesic_DistanceCot(indexMin, activeSet, activeSetDis);

			if (indexMinD < pointDistance[indexMin]) {
				pointDistance[indexMin] = indexMinD;
			}

			//add item into pointCloseDistanceHeap

			for (int i = 0; i < farSet.size(); i++) {
				float di = farSetDis[i];
				int ii = farSet[i];
				pointLabel[ii] = 1;
				pointDistance[ii] = di;
				pointCloseDistanceHeap.push_back(make_pair(di, ii));
				push_heap(pointCloseDistanceHeap.begin(), pointCloseDistanceHeap.end(), greater<pair<float, int>>());
			}

			if (iterationNum % iterstep == 0) {
				int countNum = iterationNum / iterstep;
				if (countNum == 1) {
					cout << countNum << ".";
				}
				else {
					cout << countNum << ",";
				}
			}
			//set the point to be alive and update the distance
			//quite conditions
			if (pointCloseDistanceHeap.size() <= 0) {
				judge = false;
			}
		}
		cout << endl;
		return pointDistance;

	}

	vector<int> PointGeodesic_Path(int sIndex, int tIndex) {

		vector<int> pathIndex;
		pathIndex.push_back(tIndex);

		if (sIndex >= pointSet.size() || sIndex < 0 || tIndex >= pointSet.size() || tIndex < 0) {
			cout << "Error! Input point index is out of range!";
			return pathIndex;
		}

		vector<float> dis_s = PointGeodesic_DisMap(sIndex);
		int temp = tIndex;
		double tempD = dis_s[tIndex];

		while (true) {
			int K = 7;
			std::vector<int> pointIdxNKNSearch(K);
			std::vector<float> pointNKNSquaredDistance(K);
			pcl::PointXYZ searchPoint;
			searchPoint.x = pointSet[temp][0];
			searchPoint.y = pointSet[temp][1];
			searchPoint.z = pointSet[temp][2];
			pkdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);

			int minIndex = -1;
			float minD = -9999;
			for (int i = 1; i < pointIdxNKNSearch.size(); i++) {
				float tempDi = dis_s[pointIdxNKNSearch[i]];
				float tempDii = (tempD - tempDi) / sqrt(pointNKNSquaredDistance[i]);
				if (tempDii > minD) {
					minD = tempDii;
					minIndex = pointIdxNKNSearch[i];
				}
			}
			if (minIndex == -1) {
				cout << "path searching error!";
				break;
			}
			else {
				pathIndex.push_back(minIndex);
				temp = minIndex;
				tempD = dis_s[minIndex];
			}

			if (tempD <= 0 || temp == sIndex) {
				break;
			}
		}

		return pathIndex;

	}

	vector<Point3f> PointGeodesic_Get_PointSet() {

		return pointSet;

	}

	vector<Point3f> PointGeodesic_Get_PointNormal() {

		return pointNormal;

	}

private:

	//compute the geodesic map for Non-continuous Point cloud
	void PointGeodesic_Init_PointPairDictionary(int pIndex) {

		//init global searching distance and pointLabel
		if (pointLabel.size() > 0 || pointDistance.size() > 0) {
			pointLabel.clear();
			pointDistance.clear();
		}

		if (PointPairDictionary.size() > 0) {
			PointPairDictionary.clear();
		}

		pIndexGlobal = pIndex;
		pointLabel.resize(pointSet.size(), 2);//0:Alive; 1:Close; 2:Far 
		PointPairDictionary.resize(pointSet.size());
		pointDistance.resize(pointSet.size(), 99999);
		if (pIndex >= pointLabel.size() || pIndex < 0) {
			cout << "label input error, out of vector";
			return;
		}

		int iterationNum = pointSet.size();
		vector<vector<int>> subset;
		//Achieve differnet continuous point cloud set
		int stepNum = 1;
		while (true) {
			cout << "step:" << stepNum << endl;
			stepNum++;
			vector<int> subset_i = PointGeodesic_DisMap_Unit(pIndex, iterationNum);
			subset.push_back(subset_i);
			iterationNum = iterationNum - subset_i.size();
			if (iterationNum <= 0) {
				break;
			}
			else {
				bool subJudge = true;
				for (int i = 0; i < pointLabel.size(); i++) {
					if (pointLabel[i] == 2) {
						pIndex = i;
						subJudge = false;
					}
				}
				if (subJudge) {
					break;//all points are assigned value.				
				}
			}
			cout << endl;
		}
		int maxSub = 0;
		int maxSuNum = subset[0].size();
		for (int i = 0; i < subset.size(); i++) {
			if (subset[i].size() > maxSuNum) {
				maxSuNum = subset[i].size();
				maxSub = i;
			}
		}
		//assigh class index for all point
		vector<int> pointClass(pointSet.size(), -1);

		//update subset, delete small set of point cloud part;
		vector<vector<int>> subset_new;
		for (int i = 0; i < subset.size(); i++) {
			if (i == maxSub) {
				maxSub = subset_new.size();
			}
			if (subset[i].size() > 5) {
				int classIndex = subset_new.size();
				subset_new.push_back(subset[i]);
				for (int j = 0; j < subset[i].size(); j++) {
					int indexij = subset[i][j];
					if (indexij < 0 || indexij >= pointClass.size()) {
						cout << "Point Index Error" << endl;
					}
					else {
						pointClass[indexij] = classIndex;
					}
				}
			}
		}
		subset.clear();
		subset = subset_new;
		subset_new.clear();

		/*
		//create additional neighbor searching dictionary, record the connection path to differnt parts.
		vector<vector<int>> classPath(subset.size());
		for (int i = 0; i < subset.size(); i++) {
			cout << "subset:"<< i <<endl;
			if (i == maxSub) {
				continue;
			}
			else {
				vector<pPairDic> result_i = PointGeodesic_ReturnPairDictionary(subset[i], pointClass);
				//add the pPairDic list into the PointPairDictionary
				for (int j = 0; j < result_i.size(); j++) {
					pPairDic ppd_ij = result_i[j];
					int ppd_ij_source = ppd_ij.source;
					int ppd_ij_target = ppd_ij.target;
					PointPairDictionary[ppd_ij_source].push_back(ppd_ij);
					PointPairDictionary[ppd_ij_target].push_back(ppd_ij);
					int classp1 = pointClass[ppd_ij_source];
					int classp2 = pointClass[ppd_ij_target];
					if (PointGeodesic_JudgeVector(classp2, classPath[classp1])==-1) {
						classPath[classp1].push_back(classp2);
					}
					if (PointGeodesic_JudgeVector(classp1, classPath[classp2])==-1) {
						classPath[classp2].push_back(classp1);
					}
				}
			}
			cout << endl;
		}
		*/

		if (subset.size() == 1) {
			cout << "this point cloud is a manfold point cloud!" << endl;
		}
		else {
			//Iterative clustering subset
			int clusteingIndex = 1;
			while (true) {
				cout << "Subset Clustering Iter:" << clusteingIndex << endl;
				clusteingIndex++;
				vector<vector<int>> subsetNew_Iter = PointGeodesic_Iterative_Compute_Subset(subset);
				subset.clear();
				subset = subsetNew_Iter;
				cout << "subset number:" << subset.size() << endl;
				if (subset.size() == 1) {
					break;
				}
			}
		}

		//update point cloud
		PointGeodesic_Update_Subset(subset);

		//collect different parts which is not connect to the main part

	}

	//update subset with different parts of point cloud, while adding new PointPairDictionary
	vector<vector<int>> PointGeodesic_Iterative_Compute_Subset(vector<vector<int>> subSet) {

		int maxIndex = -1;
		int maxNum = -1;
		vector<int> classCollect(subSet.size());//record the partial point clouds classification
		vector<int> pointClass(pointSet.size(), -1); // record the point classification
		for (int i = 0; i < subSet.size(); i++) {
			int subSet_i_sum = subSet[i].size();
			if (subSet_i_sum > maxNum) {
				maxIndex = i;
				maxNum = subSet[i].size();
			}
			classCollect[i] = i;
			for (int j = 0; j < subSet[i].size(); j++) {
				pointClass[subSet[i][j]] = i;
			}
		}

		vector<vector<int>> classRecord(subSet.size());

		for (int i = 0; i < subSet.size(); i++) {
			cout << "subSet:" << i << endl;
			if (i == maxIndex) {
				continue;
			}
			else {
				vector<pPairDic> result_i = PointGeodesic_ReturnPairDictionary(subSet[i], pointClass);
				//add the pPairDic list into the PointPairDictionary
				for (int j = 0; j < result_i.size(); j++) {
					pPairDic ppd_ij = result_i[j];
					int ppd_ij_source = ppd_ij.source;
					int ppd_ij_target = ppd_ij.target;
					int classp1 = pointClass[ppd_ij_source];
					int classp2 = pointClass[ppd_ij_target];
					if (classp1 == -1 || classp2 == -1) {
						continue;
					}
					PointPairDictionary[ppd_ij_source].push_back(ppd_ij);
					PointPairDictionary[ppd_ij_target].push_back(ppd_ij);

					if (PointGeodesic_JudgeVector(classp2, classRecord[classp1]) == -1) {
						classRecord[classp1].push_back(classp2);
					}
					if (PointGeodesic_JudgeVector(classp1, classRecord[classp2]) == -1) {
						classRecord[classp2].push_back(classp1);
					}
				}
			}
			cout << endl;
		}

		//check classCollect

		vector<bool> classRecordProcess(classRecord.size(), false);

		for (int i = 0; i < classRecord.size(); i++) {

			if (classRecordProcess[i]) {
				for (int j = 0; j < classRecord[i].size(); j++) {
					classRecordProcess[classRecord[i][j]] = true;
					classCollect[classRecord[i][j]] = classCollect[i];
				}
			}
			else {
				int index_Max_i = classRecord[i][0];
				int index_Sum = 0;
				for (int j = 0; j < classRecord[i].size(); j++) {
					int subSetSum_j = subSet[classRecord[i][j]].size();
					if (subSetSum_j > index_Sum) {
						index_Max_i = classRecord[i][j];
						index_Sum = subSetSum_j;
					}
				}
				for (int j = 0; j < classRecord[i].size(); j++) {
					classRecordProcess[classRecord[i][j]] = true;
					classCollect[classRecord[i][j]] = classCollect[index_Max_i];
				}
			}
		}

		//update subset
		vector<int> subSetNewIndex;
		for (int i = 0; i < classCollect.size(); i++) {
			int classCollect_i = classCollect[i];
			if (PointGeodesic_JudgeVector(classCollect_i, subSetNewIndex) == -1) {
				subSetNewIndex.push_back(classCollect_i);
			}
		}

		vector<vector<int>> subSetNew(subSetNewIndex.size());
		for (int i = 0; i < classCollect.size(); i++) {
			int classCollect_i = classCollect[i];
			int subSetNewIndex_i = PointGeodesic_JudgeVector(classCollect_i, subSetNewIndex);
			if (subSetNewIndex_i == -1) {
				cout << "Error, class searching wrong!" << endl;
			}
			else {
				subSetNew[subSetNewIndex_i].insert(subSetNew[subSetNewIndex_i].end(), subSet[i].begin(), subSet[i].end());
			}
		}

		return subSetNew;
	}

	//update new point cloud (remove outliers)
	void PointGeodesic_Update_Subset(vector<vector<int>> subset) {

		int indexCount = 0;
		vector<int> subset_dictory(pointSet.size());
		vector<Point3f> pointSet_new;
		vector<Point3f> pointNormal_new;

		for (int i = 0; i < subset.size(); i++) {
			for (int j = 0; j < subset[i].size(); j++) {
				int pIndex = subset[i][j];
				subset_dictory[pIndex] = indexCount;
				pointSet_new.push_back(pointSet[pIndex]);
				pointNormal_new.push_back(pointNormal[pIndex]);
				indexCount++;
			}
		}

		vector<vector<pPairDic>> PointPairDictionaryNew(pointSet_new.size());

		for (int i = 0; i < PointPairDictionary.size(); i++) {
			if (PointPairDictionary[i].size() > 0) {
				for (int j = 0; j < PointPairDictionary[i].size(); j++) {
					int ppd_source = PointPairDictionary[i][j].source;
					int ppd_target = PointPairDictionary[i][j].target;
					PointPairDictionary[i][j].source = subset_dictory[ppd_source];
					PointPairDictionary[i][j].target = subset_dictory[ppd_target];
				}
				PointPairDictionaryNew[subset_dictory[i]] = PointPairDictionary[i];
			}
		}

		PointPairDictionary.clear();
		PointPairDictionary = PointPairDictionaryNew;
		pointSet.clear();
		pointSet = pointSet_new;
		pointNormal.clear();
		pointNormal = pointNormal_new;
		pcl::KdTreeFLANN<pcl::PointXYZ> pkdtree_new = PointGeodesic_KdTree(pointSet);//contruct kdtree;			
		pkdtree = pkdtree_new;
	}

	//compute each continuous region in a Non-continuous Point cloud
	vector<int> PointGeodesic_DisMap_Unit(int pIndex, int iterationNum_Unit) {

		vector<int> resultPointIndex; // store the current point cloud part;
		resultPointIndex.push_back(pIndex);

		pointLabel[pIndex] = 0;
		pointDistance[pIndex] = 0;
		bool judge = true;

		//init current pIndex:
		vector<pair<float, int>> pointCloseDistanceHeap;//Close Set min-heap
		int K = 7;
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);
		pcl::PointXYZ searchPoint;
		searchPoint.x = pointSet[pIndex][0];
		searchPoint.y = pointSet[pIndex][1];
		searchPoint.z = pointSet[pIndex][2];
		pkdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);

		for (int i = 1; i < pointIdxNKNSearch.size(); i++) {
			int index_i = pointIdxNKNSearch[i];
			if (pointLabel[index_i] == 0) {
				continue;
			}
			else {
				pointLabel[index_i] = 1;
				float distance_i = sqrt(pointNKNSquaredDistance[i]);
				pointDistance[index_i] = distance_i;
				pointCloseDistanceHeap.push_back(make_pair(distance_i, index_i));
			}
		}

		if (pointCloseDistanceHeap.size() == 0) {
			cout << "the current point is an isolate outlier!" << endl;
			return resultPointIndex;
		}

		//init min-heap		
		make_heap(pointCloseDistanceHeap.begin(), pointCloseDistanceHeap.end(), greater<pair<float, int>>());
		int iterationNum = iterationNum_Unit;
		int iterstep = iterationNum / 10; //for output print
		if (iterstep == 0) {
			iterstep = 1;
		}
		//Fast Mearching
		while (judge) {
			iterationNum--;
			//0:Alive; 1:Close; 2:Far 

			//achieve the min value of the min-heap
			pop_heap(pointCloseDistanceHeap.begin(), pointCloseDistanceHeap.end(), greater<pair<float, int>>());
			pair<double, int> t = pointCloseDistanceHeap[pointCloseDistanceHeap.size() - 1];
			int indexMin = t.second;

			if (pointLabel[indexMin] == 0) {
				cout << "the current point is an isolate outlier!" << endl;
				break;
			}
			else {
				resultPointIndex.push_back(indexMin);
			}

			pointLabel[indexMin] = 0;//Set the point into Alive
			pointCloseDistanceHeap.pop_back();

			//achive the new point's neighbor
			std::vector<int> pointIdxNKNSearch_t(K);
			std::vector<float> pointNKNSquaredDistance_t(K);
			pcl::PointXYZ searchPoint_t;
			searchPoint_t.x = pointSet[indexMin][0];
			searchPoint_t.y = pointSet[indexMin][1];
			searchPoint_t.z = pointSet[indexMin][2];
			pkdtree.nearestKSearch(searchPoint_t, K, pointIdxNKNSearch_t, pointNKNSquaredDistance_t);

			vector<int> activeSet;
			vector<float> activeSetDis;
			vector<int> closeSet;
			vector<float> closeSetDis;
			vector<int> farSet;
			vector<float> farSetDis;

			//update the neighbor distance 
			for (int i = 1; i < pointIdxNKNSearch_t.size(); i++) {
				int indexi_n = pointIdxNKNSearch_t[i];
				float dis_t = sqrt(pointNKNSquaredDistance_t[i]);

				if (pointLabel[indexi_n] == 1) {
					closeSet.push_back(indexi_n);
					float closeSetDis_i = PointGeodesic_DistanceSinglePoint(indexi_n);
					if (!(closeSetDis_i <= 99999 && closeSetDis_i >= 0)) {
						cout << "Hello!" << endl;
						//pointDistance[indexi_n] = 0;
					}

					if (closeSetDis_i < pointDistance[indexi_n]) {
						pointDistance[indexi_n] = closeSetDis_i;
					}
					else {
						closeSetDis_i = pointDistance[indexi_n];
					}
					closeSetDis.push_back(closeSetDis_i);
				}
				else if (pointLabel[indexi_n] == 0) {
					activeSet.push_back(indexi_n);
					activeSetDis.push_back(pointDistance[indexi_n]);
				}
				else {
					farSet.push_back(indexi_n);
					float farSetDis_i = PointGeodesic_DistanceSinglePoint(indexi_n);
					if (farSetDis_i > dis_t + pointDistance[indexMin]) {
						farSetDis_i = dis_t + pointDistance[indexMin];
					}
					farSetDis.push_back(farSetDis_i);
				}
			}

			float indexMinD = PointGeodesic_DistanceCot(indexMin, activeSet, activeSetDis);

			if (indexMinD < pointDistance[indexMin]) {
				pointDistance[indexMin] = indexMinD;
			}

			//add item into pointCloseDistanceHeap

			for (int i = 0; i < farSet.size(); i++) {
				float di = farSetDis[i];
				int ii = farSet[i];
				pointLabel[ii] = 1;
				pointDistance[ii] = di;
				pointCloseDistanceHeap.push_back(make_pair(di, ii));
				push_heap(pointCloseDistanceHeap.begin(), pointCloseDistanceHeap.end(), greater<pair<double, int>>());
			}

			if (iterationNum % iterstep == 0) {
				int countNum = iterationNum / iterstep;
				if (countNum == 1) {
					cout << countNum << ".";
				}
				else {
					cout << countNum << ",";
				}
			}
			//set the point to be alive and update the distance
			//quite conditions
			if (pointCloseDistanceHeap.size() <= 0) {
				judge = false;
			}
		}

		return resultPointIndex;

	}

	vector<pPairDic> PointGeodesic_ReturnPairDictionary(vector<int> subset, vector<int> pointClass) {

		int K = subset.size() + 1;
		int classIndex = pointClass[subset[0]];

		vector<pPairDic> detectPairDic;
		vector<pPairDic> result;

		for (int i = 0; i < subset.size(); i++) {

			if ((subset.size() - i) % 100 == 0) {
				cout << (subset.size() - i) / 100 << ",";
			}

			int pTemIndex = subset[i];
			std::vector<int> pointIdxNKNSearch(K);
			std::vector<float> pointNKNSquaredDistance(K);
			pcl::PointXYZ searchPoint;
			searchPoint.x = pointSet[pTemIndex][0];
			searchPoint.y = pointSet[pTemIndex][1];
			searchPoint.z = pointSet[pTemIndex][2];
			pkdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);

			for (int j = 1; j < pointIdxNKNSearch.size(); j++) {
				int pTemIndex_j = pointIdxNKNSearch[j];
				if (pointClass[pTemIndex_j] != classIndex && pointClass[pTemIndex_j] != -1) {
					pPairDic ppd_ij;
					ppd_ij.source = pTemIndex_j;
					ppd_ij.target = pTemIndex;
					ppd_ij.distance = sqrt(pointNKNSquaredDistance[j]);
					detectPairDic.push_back(ppd_ij);
					break;
				}
			}
		}

		vector<bool> processedJudge(detectPairDic.size(), false);

		for (int i = 0; i < detectPairDic.size(); i++) {

			if (processedJudge[i]) {
				continue;
			}
			else {
				processedJudge[i] = true;
			}

			pPairDic ptem_i;
			ptem_i.source = detectPairDic[i].source;
			ptem_i.target = detectPairDic[i].target;
			ptem_i.distance = detectPairDic[i].distance;

			for (int j = i + 1; j < detectPairDic.size(); j++) {
				int detectPairDic_ij = detectPairDic[j].source;
				if (ptem_i.source == detectPairDic_ij) {
					processedJudge[j] = true;
					if (detectPairDic[j].distance < ptem_i.distance) {
						ptem_i.target = detectPairDic[j].target;
						ptem_i.distance = detectPairDic[j].distance;
					}
				}
			}
			result.push_back(ptem_i);
		}
		return result;
	}

	pcl::KdTreeFLANN<pcl::PointXYZ> PointGeodesic_KdTree(vector<Point3f> seedpoints) {

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

	vector<float> PointGeodesic_pointUpdate(vector<float> point_i) {

		float h;
		int iter = 10;

#pragma region Achieve Neibor
		//Achieve Neibor
		int K = 13;
		//int K = br.pointNumEsti + 1;
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);
		std::vector<int> pointNeior;
		//double r_i = pointNKNSquaredDistance[pointNKNSquaredDistance.size() - 1];
		pcl::PointXYZ searchPoint;
		searchPoint.x = point_i[0];
		searchPoint.y = point_i[1];
		searchPoint.z = point_i[2];
		pkdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);
		if (pointNKNSquaredDistance[0] == 0) {
			point_i.push_back(pointNormal[pointIdxNKNSearch[0]][0]);
			point_i.push_back(pointNormal[pointIdxNKNSearch[0]][1]);
			point_i.push_back(pointNormal[pointIdxNKNSearch[0]][2]);
			return point_i;
		}
		else {
			pointNeior.insert(pointNeior.end(), pointIdxNKNSearch.begin(), pointIdxNKNSearch.end());
		}

		//compute the h for region control in point updating
		float dis_h = sqrt(pointNKNSquaredDistance[K - 1]);
		float dis_h2 = sqrt(pointNKNSquaredDistance[K - 2]);
		h = dis_h + (dis_h - dis_h2);

		vector<Point3f> pointNormal_i(pointNeior.size());
		for (int i = 0; i < pointIdxNKNSearch.size(); i++) {
			pointNormal_i[i] = pointNormal[pointIdxNKNSearch[i]];
		}
		//regular normal
		Point3f normalUniform = pointNormal_i[0];
		for (int i = 0; i < pointNormal_i.size(); i++) {
			Point3f normalNeibor_i = pointNormal_i[i];
			Point3f normalNeibor_i2;
			normalNeibor_i2[0] = -normalNeibor_i[0];
			normalNeibor_i2[1] = -normalNeibor_i[1];
			normalNeibor_i2[2] = -normalNeibor_i[2];
			//manage different directions
			float a1 = normalUniform[0] * normalNeibor_i[0] +
				normalUniform[1] * normalNeibor_i[1] + normalUniform[2] * normalNeibor_i[2];
			float a2 = normalUniform[0] * normalNeibor_i2[0] +
				normalUniform[1] * normalNeibor_i2[1] + normalUniform[2] * normalNeibor_i2[2];
			float cosa1 = acos(a1);
			float cosa2 = acos(a2);
			if (cosa1 > cosa2) {
				pointNormal_i[i][0] = normalNeibor_i2[0];
				pointNormal_i[i][1] = normalNeibor_i2[1];
				pointNormal_i[i][2] = normalNeibor_i2[2];
			}
		}

#pragma endregion

#pragma region MLS error
		//++++++++++++++++++++interater start+++++++++++++++++++++++++++
		double errorExist = 0.0001;
		vector<float> px;
		px.insert(px.end(), point_i.begin(), point_i.end());
		//vector<int> p_neibor = br.pointNeibor[i];			
		//regularNoraml

		vector<float> px_store(3);
		vector<float> nx_store(3);
		vector<float> ax;//new point position
		ax.push_back(0);
		ax.push_back(0);
		ax.push_back(0);
		vector<float> nx;//new point normal
		nx.push_back(0);
		nx.push_back(0);
		nx.push_back(0);
		float errorEndTem;//record new 
		float errorStore = 9999;
		float weight;
		while (iter) {
			//vector<double> ax = simMeasurement_cop_a(px, pointNeiborRegualrNum);
			//vector<double> nx = simMeasurement_cop_n(px, pointNeiborRegualrNum, pointNeiborNormalRegualrNum);
			float fenmu = 0;
			for (int j = 0; j < pointNeior.size(); j++) {
				float dis_i = sqrt((pointSet[pointNeior[j]][0] - px[0]) * (pointSet[pointNeior[j]][0] - px[0]) +
					(pointSet[pointNeior[j]][1] - px[1]) * (pointSet[pointNeior[j]][1] - px[1]) +
					(pointSet[pointNeior[j]][2] - px[2]) * (pointSet[pointNeior[j]][2] - px[2]));
				if (dis_i == 0) {
					continue;
				}
				float eData = -((dis_i / h) * (dis_i / h));
				eData = exp(eData);
				ax[0] = ax[0] + pointSet[pointNeior[j]][0] * eData;
				ax[1] = ax[1] + pointSet[pointNeior[j]][1] * eData;
				ax[2] = ax[2] + pointSet[pointNeior[j]][2] * eData;
				nx[0] = nx[0] + pointNormal_i[j][0] * eData;
				nx[1] = nx[1] + pointNormal_i[j][1] * eData;
				nx[2] = nx[2] + pointNormal_i[j][2] * eData;
				fenmu = fenmu + eData;
			}
			if (fenmu != 0) {
				ax[0] = ax[0] / fenmu;
				ax[1] = ax[1] / fenmu;
				ax[2] = ax[2] / fenmu;
				nx[0] = nx[0] / fenmu;
				nx[1] = nx[1] / fenmu;
				nx[2] = nx[2] / fenmu;
			}
			fenmu = 0;

			//5.3 Set x' = x - n(x')T(a(x')-x)n(x'), weight = n(x')T(a(x')-x)
			weight = nx[0] * (point_i[0] - ax[0]) +
				nx[1] * (point_i[1] - ax[1]) + nx[2] * (point_i[2] - ax[2]);
			px_store[0] = px[0];
			px_store[1] = px[1];
			px_store[2] = px[2];
			nx_store[0] = nx[0];
			nx_store[1] = nx[1];
			nx_store[2] = nx[2];
			px[0] = point_i[0] - weight * nx[0];
			px[1] = point_i[1] - weight * nx[1];
			px[2] = point_i[2] - weight * nx[2];
			//5.4 ||n(x')T(a(x')-x)n(x')||>errorEndTem
			ax[0] = 0;
			ax[1] = 0;
			ax[2] = 0;
			nx[0] = 0;
			nx[1] = 0;
			nx[2] = 0;

			//vector<double> axnew = simMeasurement_cop_a(px, pointNeiborRegualrNum);
			//vector<double> nxnew = simMeasurement_cop_n(px, pointNeiborRegualrNum, pointNeiborNormalRegualrNum);

			for (int j = 0; j < pointNeior.size(); j++) {
				float dis_i = sqrt((pointSet[pointNeior[j]][0] - px[0]) * (pointSet[pointNeior[j]][0] - px[0]) +
					(pointSet[pointNeior[j]][1] - px[1]) * (pointSet[pointNeior[j]][1] - px[1]) +
					(pointSet[pointNeior[j]][2] - px[2]) * (pointSet[pointNeior[j]][2] - px[2]));
				if (dis_i == 0) {
					continue;
				}
				float eData = -((dis_i / h) * (dis_i / h));
				eData = exp(eData);
				ax[0] = ax[0] + pointSet[pointNeior[j]][0] * eData;
				ax[1] = ax[1] + pointSet[pointNeior[j]][1] * eData;
				ax[2] = ax[2] + pointSet[pointNeior[j]][2] * eData;
				nx[0] = nx[0] + pointNormal_i[j][0] * eData;
				nx[1] = nx[1] + pointNormal_i[j][1] * eData;
				nx[2] = nx[2] + pointNormal_i[j][2] * eData;
				fenmu = fenmu + eData;
			}
			if (fenmu != 0) {
				ax[0] = ax[0] / fenmu;
				ax[1] = ax[1] / fenmu;
				ax[2] = ax[2] / fenmu;
				nx[0] = nx[0] / fenmu;
				nx[1] = nx[1] / fenmu;
				nx[2] = nx[2] / fenmu;
			}

			weight = nx[0] * (point_i[0] - ax[0]) +
				nx[1] * (point_i[1] - ax[1]) + nx[2] * (point_i[2] - ax[2]);

			ax[0] = 0;
			ax[1] = 0;
			ax[2] = 0;
			nx[0] = 0;
			nx[1] = 0;
			nx[2] = 0;

			errorEndTem = abs(weight);
			if (errorEndTem < errorStore) {
				errorStore = errorEndTem;
			}
			else {
				px[0] = px_store[0];
				px[1] = px_store[1];
				px[2] = px_store[2];
				break;
			}
			if (errorEndTem < errorExist) {//|| errorEndTem > errorstore
				break;
			}
			iter--;
		}
#pragma endregion 

		vector<float> finalResult(6);
		finalResult[0] = px[0];
		finalResult[1] = px[1];
		finalResult[2] = px[2];
		finalResult[3] = nx_store[0];
		finalResult[4] = nx_store[1];
		finalResult[5] = nx_store[2];
		return finalResult;
	}

	float PointGeodesic_DistanceSinglePoint(int pIndex) {

		int K = 7;
		std::vector<int> pointIdxNKNSearch_t(K);
		std::vector<float> pointNKNSquaredDistance_t(K);
		pcl::PointXYZ searchPoint_t;
		searchPoint_t.x = pointSet[pIndex][0];
		searchPoint_t.y = pointSet[pIndex][1];
		searchPoint_t.z = pointSet[pIndex][2];
		pkdtree.nearestKSearch(searchPoint_t, K, pointIdxNKNSearch_t, pointNKNSquaredDistance_t);

		vector<int> activeSetSP;
		vector<float> activeSetDisSP;

		for (int i = 1; i < pointIdxNKNSearch_t.size(); i++) {
			int indexi = pointIdxNKNSearch_t[i];
			if (pointLabel[indexi] == 0) {
				activeSetSP.push_back(indexi);
				activeSetDisSP.push_back(pointDistance[indexi]);
			}
		}

		float pIndexDis = PointGeodesic_DistanceCot(pIndex, activeSetSP, activeSetDisSP);
		return pIndexDis;

	}

	float PointGeodesic_DistanceCot(int pIndex, vector<int> activeSet, vector<float> activeSetDis) {

		int x1Index;
		int x2Index;
		float a1;
		float a2;

		if (activeSet.size() <= 1) {
			//if (activeSet.size() == 1) {
				//return activeSetDis[0];			
			//}
			//else {
			return 9999;
			//}			
		}
		else if (activeSet.size() == 2) {
			x1Index = activeSet[0];
			x2Index = activeSet[1];
			a1 = activeSetDis[0];
			a2 = activeSetDis[1];

		}
		else {
			int waitSet;
			x1Index = activeSet[0];
			a1 = activeSetDis[0];
			for (int i = 0; i < activeSet.size(); i++) {
				if (activeSetDis[i] < a1) {
					a1 = activeSetDis[i];
					x1Index = activeSet[i];
				}
			}

			if (x1Index == activeSet[0]) {
				waitSet = 1;
			}
			else {
				waitSet = 0;
			}

			int K = 13;
			//int K = br.pointNumEsti + 1;
			std::vector<int> pointIdxNKNSearch(K);
			std::vector<float> pointNKNSquaredDistance(K);
			std::vector<int> pointNeior;
			//double r_i = pointNKNSquaredDistance[pointNKNSquaredDistance.size() - 1];
			pcl::PointXYZ searchPoint;
			searchPoint.x = pointSet[x1Index][0];
			searchPoint.y = pointSet[x1Index][1];
			searchPoint.z = pointSet[x1Index][2];
			pkdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);

			int indexSearch = -1;
			for (int i = 1; i < pointIdxNKNSearch.size(); i++) {
				if (pointIdxNKNSearch[i] != pIndex && pointIdxNKNSearch[i] != x1Index) {
					indexSearch = PointGeodesic_JudgeVector(pointIdxNKNSearch[i], activeSet);
					if (indexSearch >= 0) {
						break;
					}
				}
			}

			if (indexSearch = -1) {
				x2Index = activeSet[waitSet];
				a2 = activeSetDis[waitSet];
			}
			else {
				x2Index = activeSet[indexSearch];
				a2 = activeSetDis[indexSearch];
			}

		}

		if (a1 == 0) {
			return PointGeodesic_DistancePoint(pIndex, x1Index);

		}

		if (a2 == 0) {
			return PointGeodesic_DistancePoint(pIndex, x2Index);
		}

		//if (a1 == a2) {
			//return a1;
		//}

		//x1: (0,0)
		//x2: (0,d(x1,x2))
		//x3: (d(x1,x3)*cosA, d(x1,x3)*sinA)
		//x4: ()
		//x3:
		float bt = PointGeodesic_DistancePoint(x1Index, pIndex);
		float ct = PointGeodesic_DistancePoint(x2Index, pIndex);
		float at = PointGeodesic_DistancePoint(x2Index, x1Index);;
		float cosA = (bt * bt + ct * ct - at * at) / (2 * bt * ct);
		if (cosA > 1) {
			cosA = 1;
		}
		if (cosA < -1) {
			cosA = -1;
		}
		float sinA = sin(acos(cosA));
		float x3 = abs(PointGeodesic_DistancePoint(x1Index, pIndex) * cosA);
		float y3 = abs(PointGeodesic_DistancePoint(x1Index, pIndex) * sinA);
		//x4:
		float x4 = -(a2 * a2 - a1 * a1 - at * at) / (2 * at);
		float y4_f = abs(a1 * a1 - x4 * x4);
		float y4;
		if (y4_f < 0.000001) {
			y4 = 0;
		}
		else {
			y4 = -sqrt(y4_f);
		}
		float result = sqrt((x3 - x4) * (x3 - x4) + (y3 - y4) * (y3 - y4));
		return result;

	}

	float PointGeodesic_DistancePoint(int p1, int p2) {

		Point3f pIndex1 = pointSet[p1];
		Point3f pIndex2 = pointSet[p2];
		float pDis = sqrt((pIndex1[0] - pIndex2[0]) * (pIndex1[0] - pIndex2[0]) +
			(pIndex1[1] - pIndex2[1]) * (pIndex1[1] - pIndex2[1]) +
			(pIndex1[2] - pIndex2[2]) * (pIndex1[2] - pIndex2[2]));

		return pDis;

	}

	int PointGeodesic_JudgeVector(int p1, vector<int> t) {

		for (int i = 0; i < t.size(); i++) {

			if (t[i] == p1) {
				return i;

			}

		}
		return -1;

	}

};