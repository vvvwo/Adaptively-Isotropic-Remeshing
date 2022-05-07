/********************************************************************************
*
*                           SHFRemshing Main Function
*
*                               by Dr. Chenlei Lv
*
*                                  2021.08.02
*
*********************************************************************************/
#pragma once
#include "Remeshing/SHFRemesh.cpp"
using namespace std;

int main(int argc, char** argv)
{
    cout << "start!" << endl;
    char* pfile;     
    char* preprocessing;
    char* isotropicJudge;
    char* meshscaleInput;
    int curvatureSensative;
    if (argc != 5) {
        cout << "Parameters input error!" << endl;      
        //string filePath = "Data/Mesh/hand.obj";
        //string filePath = "Data/Mesh/bowl.obj";
        //SHFRemeshing shr;
        //shr.SHFRemeshing_init(filePath);
        //shr.SHFRemeshing_Start(2);
        //shr.SHFRemeshing_Start_Pro(2, 0.2);
    }
    else {
        pfile = argv[1];        
        preprocessing = argv[2];
        isotropicJudge = argv[3];
        meshscaleInput = argv[4];

        //pfile = "E:/chen_database/RemshOriginal/stand/Bunny.obj";
        //pfile = "D:/chen/research/project/Active_Anisotropic Remesh/data/_Inital/car_0089.obj";
        //pfile = "E:/download/model/objList/experimental/Lu Yu-obj_iso.obj";
        //pfile = "E:/chen_database/RemshOriginal/shrec/hand.obj";
        //pfile = "E:/chen_database/_Remesh/TestData/SHREC/S_tyrannosaurusT69.obj";
        //pfile = "E:/chen_database/_Remesh/TestData/experimental/Garuda_and_Vishnu_iso.obj";
        //preprocessing = "0"; //0: do not process; 1: process
        //isotropicJudge = "2"; //1: isotropic 2: An-isotropic
        //meshscaleInput = "1.2"; //
        curvatureSensative = 1;

        int preprocessingindex = atoi(preprocessing);
        int isotropicJudgeindex = atoi(isotropicJudge);
        double meshscale = atof(meshscaleInput);
        string filePath1(pfile);
        SHFRemeshing shr;
        shr.SHFRemeshing_init(filePath1);
        //shr.SHFRemeshing_Start(2);
        if (preprocessingindex == 1) {
            shr.SHFRemeshing_Start_Pro(isotropicJudgeindex, curvatureSensative, meshscale);           
        }
        else {
            shr.SHFRemeshing_Start(isotropicJudgeindex, curvatureSensative, meshscale);
        }
    }    
    
}




