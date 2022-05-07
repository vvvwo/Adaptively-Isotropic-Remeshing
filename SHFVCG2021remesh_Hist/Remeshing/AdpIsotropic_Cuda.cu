#include "cuda_runtime.h"  
#include "device_launch_parameters.h"  
#include <vector>

__global__ void AnIsotropic_Unit_Cuda_Block(
	float** ps_2d_Cu,
	float* pt_para_Cu,
	int* XYZNumber_Cu,
	int* pt_VoxelJudge_Cu,
	int voxelSize,
	int length,//ps_2d_Cu length for a center
	int centerNum,//source center number
	float rate,
	float* resultMatching_Cu
) {

}

extern "C" void AnIsotropic_Unit_Cuda(
	std::vector<std::vector<double>> points,//source point cloud with different centers	
	std::vector<std::vector<int>> pointNeighbor,//parameter of target point cloud voxel structure: minXYZ,maxXYZ,unitSize
	std::vector<double> cu_ave,//parameter of target point cloud voxel structure: XYZNumber
	std::vector<int> points_Keep
) {



}