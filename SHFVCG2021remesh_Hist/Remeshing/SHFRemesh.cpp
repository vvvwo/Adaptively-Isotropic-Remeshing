/********************************************************************************
*
*                                SHFRemesh
*
*                             by Dr. Chenlei Lv
*
*                                2021.10.21
*
*********************************************************************************/
#pragma once

#include "TriangleCheck.cpp"
#include "TriangleRemoveRepeat.cpp"
#include "SubRemeshing.cpp"
#include "AdpIsotropic.cpp"
#include "AdpIsotropic_OM.cpp"
#include "Isotropic.cpp"
#include "../Load_Resampling.hpp"
#include "../Normal_Generate.hpp"
#include "../IntrinsicTransfer.hpp"


class SHFRemeshing {

private:

	string fileNameInput;//store basic path
	string fileNameObj;//store basic path
	double h;
	double weightPointUpdate = 0.6;//update weight for seedpoints
	
	vector<vector<double>> pointSet;
	vector<vector<int>> faceSet;
	double meshScale;//remesh, scaling of average edges
	int aptCuvature;

	vector<vector<double>> pointSet_Original;
	vector<vector<double>> pointSet_Original_Normal;
	pcl::KdTreeFLANN<pcl::PointXYZ> pointSet_Original_KD;

public:

	void SHFRemeshing_init(string fileName) {
		
		char* p = new char[strlen(fileName.c_str()) + 1];
		strcpy(p, fileName.c_str());
		Load_Resampling lr;
		lr.Load_Resampling_init(p);
		vector<Point3f> point3f = lr.pointOriginal;		
		pointSet.resize(point3f.size());
		faceSet = lr.faceOriginal;
		for (int i = 0; i < point3f.size(); i++) {
			pointSet[i].push_back((double)point3f[i][0]);
			pointSet[i].push_back((double)point3f[i][1]);
			pointSet[i].push_back((double)point3f[i][2]);		
		}
		pointSet_Original = pointSet;

		//init KDTree		
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->width = pointSet_Original.size();
		cloud->height = 1;
		cloud->points.resize(cloud->width * cloud->height);
		// fills a PointCloud with random data
		for (int i = 0; i < pointSet_Original.size(); i++)
		{
			pcl::PointXYZ pxyz;
			cloud->points[i].x = pointSet_Original[i][0];
			cloud->points[i].y = pointSet_Original[i][1];
			cloud->points[i].z = pointSet_Original[i][2];

		}
		pointSet_Original_KD.setInputCloud(cloud);

		//init normal
		Normal_Generate ng;
		vector<Point3f> normalVector = ng.estimateNormal_PCL_MP(point3f, 8);
		for (int i = 0; i < normalVector.size();i++) {
			vector<double> pi(3);
			pi[0] = normalVector[i][0];
			pi[1] = normalVector[i][1];
			pi[2] = normalVector[i][2];
			pointSet_Original_Normal.push_back(pi);		
		}

		int index = fileName.find_last_of(".");
		fileNameInput = fileName.substr(0, index);
		int indexName = fileName.find_last_of("/");
		fileNameObj = fileName.substr(indexName+1, index-indexName-1);
	
	}

	void SHFRemeshing_Start(int classification, int curvaturesensative, double meshScaleInput) {//the classification determine the isotropic strategy
		//classificationIndex 1: isotropic; 2: Curvature; 3: Edge.
		
		meshScale = meshScaleInput;
		aptCuvature = curvaturesensative;
		cout << "SHFRemeshing Start:" << endl;
		//Pre-processing
		//SHFRemeshing_TriangleRemoveRepeat();
		//SHFRemeshing_TriangleCheck();
		//SHFRemeshing_SubRemeshing();

		if (classification == 1) {
			SHFRemeshing_Isotropic();		
		}		
		else if (classification == 2) {
			SHFRemeshing_AnIsotropic();		
		}
		else {
			cout << "Parameter Error!" << endl;		
		}		

	}

	void SHFRemeshing_Start_Pro(int classification, int curvaturesensative, double meshScaleInput) {//the classification determine the isotropic strategy
		//classificationIndex 1: isotropic; 2: Curvature; 3: Edge.

		meshScale = meshScaleInput;
		aptCuvature = curvaturesensative;
		cout << "SHFRemeshing Start:" << endl;
		//Pre-processing
		SHFRemeshing_TriangleRemoveRepeat();
		SHFRemeshing_TriangleCheck();
		SHFRemeshing_SubRemeshing();

		if (classification == 1) {
			SHFRemeshing_Isotropic();
		}
		else if (classification == 2) {
			SHFRemeshing_AnIsotropic();
		}
		else {
			cout << "Parameter Error!" << endl;
		}

	}

private:

	void SHFRemeshing_TriangleRemoveRepeat() {

		Triangluar_RemoveRepeat tr;
		tr.TR_Init(pointSet, faceSet);
		tr.TR_Start(fileNameObj);
		pointSet.clear();
		faceSet.clear();
		pointSet = tr.TR_Get_PointSet();
		faceSet = tr.TR_Get_FaceSet();
	
	}

	void SHFRemeshing_TriangleCheck() {

		Triangluar_Check tc;
		tc.TC_init(pointSet, faceSet);
		tc.TC_start(fileNameObj, false);
		pointSet.clear();
		faceSet.clear();
		pointSet = tc.TC_Get_pointSet();
		faceSet = tc.TC_Get_faceSet();

	}

	void SHFRemeshing_SubRemeshing() {

		Sub_Remeshing sr;
		sr.Sub_Remeshing_init(pointSet, faceSet, 1.5);
		sr.Sub_Remeshing_Start(fileNameObj);
		pointSet.clear();
		faceSet.clear();
		pointSet = sr.Sub_Get_PointSet();
		faceSet = sr.Sub_Get_FaceSet();

	}

	void SHFRemeshing_Isotropic() {

		Isotropic mo;
		mo.MeshOptimization_init(fileNameObj, pointSet, faceSet);
		int t0 = clock();		
		mo.MeshOptimization_Start(meshScale, 5, false);
		int t1 = clock();
		cout << ",finished! time:" << float(t1 - t0) / CLOCKS_PER_SEC << "s" << endl;
		pointSet.clear();
		pointSet = mo.MeshOptimization_GetPoints();
		faceSet.clear();
		faceSet = mo.MeshOptimization_GetFaces();
		MeshOptimization_SaveOBJ_Face(pointSet, faceSet);
	
	}

	void SHFRemeshing_AnIsotropic() {

		//The second parameter: AnIsotropic Classification
	    //1:cu_ave[i] = {1.2, 1.1, 1.0, 0.9. 0.8}
	    //2:cu_ave[i] = {3.0, 2.5, 2, 1.5. 1.2, 1, 0.8, 0.7}
		//Isotropic mo1;
		//mo1.MeshOptimization_init(fileNameObj, pointSet, faceSet);
		//mo1.MeshOptimization_Start(meshScale, 5, true);
		//pointSet.clear();
		//pointSet = mo1.MeshOptimization_GetPoints();
		//faceSet.clear();
		//faceSet = mo1.MeshOptimization_GetFaces();

		
		AnIsotropic mo;
		mo.MeshOptimization_init(fileNameObj, aptCuvature, pointSet, faceSet);//2 parameter: kind of adapter
		int t0 = clock();
		mo.MeshOptimization_Start(meshScale, 5, false);//2 parameter: iter number
		int t1 = clock();
		cout << ",finished! time:" << float(t1 - t0) / CLOCKS_PER_SEC << "s" << endl;
		pointSet.clear();
		pointSet = mo.MeshOptimization_GetPoints();
		faceSet.clear();
		faceSet = mo.MeshOptimization_GetFaces();

		//align pointSet
		cout << "Update point location" << endl;
		for (int i = 0; i < pointSet.size(); i++) {

			vector<double> pointSetNew = MeshOptimization_PointUpdate(pointSet[i]);
			pointSet[i] = pointSetNew;
		
		}

		MeshOptimization_SaveOBJ_Face(pointSet, faceSet);				
		

		/*
		AnIsotropic_OM moo;
		moo.AnIsotropic_OM_init(fileNameObj, aptCuvature, pointSet, faceSet);
		int t0 = clock();
		moo.AnIsotropic_OM_Start(meshScale, 5, false);
		int t1 = clock();
		cout << ",finished! time:" << float(t1 - t0) / CLOCKS_PER_SEC << "s" << endl;
		pointSet.clear();
		pointSet = moo.MeshOptimization_GetPoints();
		faceSet.clear();
		faceSet = moo.MeshOptimization_GetFaces();
		MeshOptimization_SaveOBJ_Face(pointSet, faceSet);
		*/

	}

	void MeshOptimization_SaveOBJ_Face(vector<vector<double>> points, vector<vector<int>> faceInfor) {

		string fin = fileNameInput + "_Remesh_"+ to_string(meshScale) +".obj";

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

	//vector<vector<double>> pointSet_Original;
	//vector<vector<double>> pointSet_Original_Normal;
	//pcl::KdTreeFLANN<pcl::PointXYZ> pointSet_Original_KD;
	vector<double> MeshOptimization_PointUpdate(vector<double> point_i) {

		int iter = 10;

#pragma region Achieve Neibor
		//Achieve Neibor
		int K = 8;
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);
		std::vector<int> pointNeior;
		//double r_i = pointNKNSquaredDistance[pointNKNSquaredDistance.size() - 1];
		pcl::PointXYZ searchPoint;
		searchPoint.x = point_i[0];
		searchPoint.y = point_i[1];
		searchPoint.z = point_i[2];
		pointSet_Original_KD.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);
		if (pointNKNSquaredDistance[0] == 0) {
			point_i.push_back(pointSet_Original_Normal[pointIdxNKNSearch[0]][0]);
			point_i.push_back(pointSet_Original_Normal[pointIdxNKNSearch[0]][1]);
			point_i.push_back(pointSet_Original_Normal[pointIdxNKNSearch[0]][2]);
			return point_i;
		}
		else {
			pointNeior.insert(pointNeior.end(), pointIdxNKNSearch.begin(), pointIdxNKNSearch.end() - 1);
		}


		//vector<int> pointNeior = br.pointNeibor[i];
		vector<vector<double>> pointNormal_i(pointNeior.size());
		pointNormal_i[0] = pointSet_Original_Normal[pointNeior[0]];
		for (int j = 1; j < pointNeior.size(); j++) {
			//vector<double> n_j = br.pointNormal[pointNeior[j]];
			pointNormal_i[j] = pointSet_Original_Normal[pointNeior[j]];
		}

#pragma endregion

#pragma region MLS error
		//++++++++++++++++++++interater start+++++++++++++++++++++++++++
		double errorExist = 0.0001;
		vector<double> px;
		px.insert(px.end(), point_i.begin(), point_i.end());
		//vector<int> p_neibor = br.pointNeibor[i];			
		//regularNoraml

		vector<double> px_store(3);
		vector<double> nx_store(3);
		vector<double> ax;//new point position
		ax.push_back(0);
		ax.push_back(0);
		ax.push_back(0);
		vector<double> nx;//new point normal
		nx.push_back(0);
		nx.push_back(0);
		nx.push_back(0);
		double errorEndTem;//record new 
		double errorStore = 9999;
		double weight;
		while (iter) {
			//vector<double> ax = simMeasurement_cop_a(px, pointNeiborRegualrNum);
			//vector<double> nx = simMeasurement_cop_n(px, pointNeiborRegualrNum, pointNeiborNormalRegualrNum);
			double fenmu = 0;
			for (int j = 0; j < pointNeior.size(); j++) {
				double dis_i = sqrt((pointSet_Original[pointNeior[j]][0] - px[0]) * (pointSet_Original[pointNeior[j]][0] - px[0]) +
					(pointSet_Original[pointNeior[j]][1] - px[1]) * (pointSet_Original[pointNeior[j]][1] - px[1]) +
					(pointSet_Original[pointNeior[j]][2] - px[2]) * (pointSet_Original[pointNeior[j]][2] - px[2]));
				if (dis_i == 0) {
					continue;
				}
				double eData = -((dis_i / h) * (dis_i / h));
				eData = exp(eData);
				ax[0] = ax[0] + pointSet_Original[pointNeior[j]][0] * eData;
				ax[1] = ax[1] + pointSet_Original[pointNeior[j]][1] * eData;
				ax[2] = ax[2] + pointSet_Original[pointNeior[j]][2] * eData;
				nx[0] = nx[0] + pointSet_Original_Normal[pointNeior[j]][0] * eData;
				nx[1] = nx[1] + pointSet_Original_Normal[pointNeior[j]][1] * eData;
				nx[2] = nx[2] + pointSet_Original_Normal[pointNeior[j]][2] * eData;
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
				double dis_i = sqrt((pointSet_Original[pointNeior[j]][0] - px[0]) * (pointSet_Original[pointNeior[j]][0] - px[0]) +
					(pointSet_Original[pointNeior[j]][1] - px[1]) * (pointSet_Original[pointNeior[j]][1] - px[1]) +
					(pointSet_Original[pointNeior[j]][2] - px[2]) * (pointSet_Original[pointNeior[j]][2] - px[2]));
				if (dis_i == 0) {
					continue;
				}
				double eData = -((dis_i / h) * (dis_i / h));
				eData = exp(eData);
				ax[0] = ax[0] + pointSet_Original[pointNeior[j]][0] * eData;
				ax[1] = ax[1] + pointSet_Original[pointNeior[j]][1] * eData;
				ax[2] = ax[2] + pointSet_Original[pointNeior[j]][2] * eData;
				nx[0] = nx[0] + pointSet_Original_Normal[pointNeior[j]][0] * eData;
				nx[1] = nx[1] + pointSet_Original_Normal[pointNeior[j]][1] * eData;
				nx[2] = nx[2] + pointSet_Original_Normal[pointNeior[j]][2] * eData;
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

		vector<double> finalResult(3);
		finalResult[0] = px[0];
		finalResult[1] = px[1];
		finalResult[2] = px[2];
		//finalResult[3] = nx_store[0];
		//finalResult[4] = nx_store[1];
		//finalResult[5] = nx_store[2];
		return finalResult;

	}

};



