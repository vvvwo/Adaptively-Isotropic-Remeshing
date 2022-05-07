/*********************************************************************
*
*                   Generate Point Normal Vectors
*
*                         by Dr. Chenlei Lv
*
*                            2021.08.06
*
**********************************************************************/
#include <vector>
#include <vcg/complex/complex.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d_omp.h>

using namespace std;
using namespace vcg;

class Normal_Generate {	
	
public:	

	vector<Point3f> estimateNormal_PCL_MP(vector<Point3f> pointsVector, int neighborNum) {

		vector<Point3f> normalPoint;

		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);

		// fills a PointCloud with random data
		for (int i = 0; i < pointsVector.size(); i++)
		{
			pcl::PointXYZ cloud_i;
			cloud_i.x = pointsVector[i][0];
			cloud_i.y = pointsVector[i][1];
			cloud_i.z = pointsVector[i][2];
			cloud->push_back(cloud_i);
		}

		// Create the normal estimation class, and pass the input dataset to it
		pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne;
		ne.setNumberOfThreads(12);
		ne.setInputCloud(cloud);

		// Create an empty kdtree representation, and pass it to the normal estimation object.
		// Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
		pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
		ne.setSearchMethod(tree);
		// Output datasets
		pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(new pcl::PointCloud<pcl::Normal>);

		// Use all neighbors in a sphere of radius 3cm
		ne.setKSearch(neighborNum);
		//ne.setRadiusSearch(0.03);

		// Compute the features
		ne.compute(*cloud_normals);

		for (int i = 0; i < cloud_normals->size(); i++) {

			Point3f normal_i;
			float dis_i = sqrt((cloud_normals->at(i).normal_x) * (cloud_normals->at(i).normal_x) +
				(cloud_normals->at(i).normal_y) * (cloud_normals->at(i).normal_y) +
				(cloud_normals->at(i).normal_z) * (cloud_normals->at(i).normal_z));
			normal_i[0] = cloud_normals->at(i).normal_x / dis_i;
			normal_i[1] = cloud_normals->at(i).normal_y / dis_i;
			normal_i[2] = cloud_normals->at(i).normal_z / dis_i;
			normalPoint.push_back(normal_i);

		}

		return normalPoint;

	}





};