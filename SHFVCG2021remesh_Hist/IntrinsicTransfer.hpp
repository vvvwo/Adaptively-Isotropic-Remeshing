/*********************************************************************************

				Intrinsic transfer pro for a non-manifold point cloud

							 Updating in 2021/06/29

							   By Dr. Chenlei Lv

			The functions includes:
			1. Using the geodesic distance field to transfer a point cloud into
			   an intrinsic one.



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
#include "geodesic.hpp"

using namespace std;
using namespace vcg;

class IntrinsicTransfer_Pro {

public:

	vector<float> geodesicDistance;
	float gMax;
	float gMin;

private:

	vector<Point3f> pointSet;
	vector<Point3f> pointSetNormal;
	vector<Point3f> pointSerIntrinsic;
	PointGeodesic_Pro pg_pro;

public:

	void IntrinsicTransfer_init(vector<Point3f> pointSet_input, vector<Point3f> pointSetNormal_input) {

		pg_pro.PointGeodesic_Pro_init(pointSet_input, pointSetNormal_input);
		pointSet = pg_pro.PointGeodesic_Get_PointSet();
		pointSetNormal = pg_pro.PointGeodesic_Get_PointNormal();

	}

	vector<Point3f> IntrinsicTransfer_IntrinsicMap() {

		IntrinsicTransfer_geoMapping();
		IntrinsicTransfer_Uniform();
		return pointSerIntrinsic;

	}

	vector<Point3f> IntrinsicTransfer_Get_PointSet() {

		return pointSet;

	}

	vector<Point3f> IntrinsicTransfer_Get_PointNormal() {

		return pointSetNormal;

	}

private:

	void IntrinsicTransfer_geoMapping() {

		//non-manifold point cloud
	//PointGeodesic_Pro pg;
	//pg.PointGeodesic_Pro_init(ppl.br.pointCloudData, ppl.br.pointNormal);
	//distanceG = pg.PointGeodesic_DisMap(0);
	//vector<vector<double>> pointSet = pg.PointGeodesic_Get_PointSet();

		gMax = -9999;
		gMin = 9999;

		int tempIndex = 0;
		vector<float> dis = pg_pro.PointGeodesic_DisMap(tempIndex);
		float maxDis = -1;
		int firstPoint = -1;
		for (int i = 0; i < dis.size(); i++) {
			if (dis[i] > maxDis) {
				maxDis = dis[i];
				firstPoint = i;
			}
		}

		vector<float> disFirst = pg_pro.PointGeodesic_DisMap(firstPoint);
		maxDis = -1;
		int secondPoint = -1;
		for (int i = 0; i < disFirst.size(); i++) {
			if (disFirst[i] > maxDis) {
				maxDis = disFirst[i];
				secondPoint = i;
			}
		}

		vector<float> disSecond = pg_pro.PointGeodesic_DisMap(secondPoint);

		maxDis = -1;
		int thirdPoint = -1;
		for (int i = 0; i < dis.size(); i++) {
			float dis3i;
			if (disFirst[i] > disSecond[i]) {
				dis3i = disSecond[i];
			}
			else {
				dis3i = disFirst[i];
			}
			if (dis3i > maxDis) {
				maxDis = dis3i;
				thirdPoint = i;
			}
		}

		vector<float> disThird = pg_pro.PointGeodesic_DisMap(thirdPoint);

		if (pointSerIntrinsic.size() > 0) {
			pointSerIntrinsic.clear();
		}

		geodesicDistance.clear();
		geodesicDistance.resize(dis.size(), 0);

		for (int i = 0; i < dis.size(); i++) {
			float xi = disFirst[i];
			float yi = disSecond[i];
			float zi = disThird[i];
			Point3f pointIndex;
			pointIndex[0] = xi;
			pointIndex[1] = yi;
			pointIndex[2] = zi;
			pointSerIntrinsic.push_back(pointIndex);

			if (xi > yi) {
				xi = yi;
			}
			if (xi > zi) {
				xi = zi;
			}

			geodesicDistance[i] = xi;
			if (xi > gMax) {
				gMax = xi;
			}
			if (xi < gMin) {
				gMin = xi;
			}
		}

	}

	void IntrinsicTransfer_Uniform() {

		float maxx = -9999, minx = 9999;
		float maxy = -9999, miny = 9999;
		float maxz = -9999, minz = 9999;

		for (int i = 0; i < pointSerIntrinsic.size(); i++) {

			float xi = pointSerIntrinsic[i][0];
			float yi = pointSerIntrinsic[i][1];
			float zi = pointSerIntrinsic[i][2];

			if (xi > maxx) {
				maxx = xi;
			}
			if (xi < minx) {
				minx = xi;
			}
			if (yi > maxy) {
				maxy = yi;
			}
			if (yi < miny) {
				miny = yi;
			}
			if (zi > maxz) {
				maxz = zi;
			}
			if (zi < minz) {
				minz = zi;
			}

		}

		float xCenter = (maxx - minx) / 2;
		float yCenter = (maxy - miny) / 2;
		float zCenter = (maxz - minz) / 2;


		for (int i = 0; i < pointSerIntrinsic.size(); i++) {

			pointSerIntrinsic[i][0] = pointSerIntrinsic[i][0] - xCenter;
			pointSerIntrinsic[i][1] = pointSerIntrinsic[i][1] - yCenter;
			pointSerIntrinsic[i][2] = pointSerIntrinsic[i][2] - zCenter;

		}

	}
};