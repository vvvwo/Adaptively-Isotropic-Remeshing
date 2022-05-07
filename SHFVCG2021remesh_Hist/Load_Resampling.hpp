/*********************************************************************
* 
*         Load and Resampling Point Cloud from input file       
* 
*                     by Dr. Chenlei Lv
*
*                          2021.08.04 
*  
**********************************************************************/

#pragma once
#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/clustering.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
    Use<MyEdge>     ::AsEdgeType,
    Use<MyFace>     ::AsFaceType> {};

class MyVertex : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::BitFlags, vertex::Mark> {};
class MyFace : public Face< MyUsedTypes, face::Mark, face::VertexRef, face::VFAdj, face::FFAdj, face::Normal3f, face::BitFlags > {};
class MyEdge : public Edge<MyUsedTypes> {};
class MyMesh : public tri::TriMesh< vector<MyVertex>, vector<MyFace>, vector<MyEdge>> {};

class Load_Resampling {

public:

    float radiusGlobal = 0;
    MyMesh m;  
    vector<Point3f> pointOriginal;
    vector<vector<int>> faceOriginal;

public:

    void Load_Resampling_init(char* pfile) {       
       
        //load mesh
        cout << "Load data...";
        //int t0 = clock();
        if (tri::io::Importer<MyMesh>::Open(m, pfile) != 0)
        {
            printf("Error reading file  %s\n", pfile);
            exit(0);
        }
        //int t1 = clock();
        Load_Resampling_Store();
        //cout << ",finished! time:" << float(t1 - t0) / CLOCKS_PER_SEC << "s" << endl;   
        cout << endl;
    
    }

    vector<Point3f> Load_Resampling_Simplify(int simplification) {
       
       // int t0 = clock();
        MyMesh subM;
        tri::MeshSampler<MyMesh> mps(subM);        

        if (faceOriginal.size() < 1) {
            
            //cout << "Estimate Possion Radius:";
            float radius = tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh>>::ComputePoissonDiskRadius(m, simplification * 1.2);
            radiusGlobal = radius;
            //cout << radius << endl;
            cout << "Pre-processing..." << endl;
            tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh>>::PoissonDiskParam pp;
            pp.bestSampleChoiceFlag = false;
            tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh>>::PoissonDiskPruning(mps, m, radius, pp);
            //int t1 = clock();
           // cout << ",finished! Sampling points:" << subM.VN() << "time:" << float(t1 - t0) / CLOCKS_PER_SEC << "s" << endl;           
        
        }
        else {
            MyMesh m_wFace;
            tri::BuildMeshFromCoordVector(m_wFace, pointOriginal);
            float radius = tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh>>::ComputePoissonDiskRadius(m_wFace, simplification * 1.2);
            radiusGlobal = radius;
            //cout << radius << endl;
            cout << "Pre-processing..." << endl;
            tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh>>::PoissonDiskParam pp;
            pp.bestSampleChoiceFlag = false;
            tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh>>::PoissonDiskPruning(mps, m_wFace, radius, pp);
           // int t1 = clock();
            //cout << ",finished! Sampling points:" << subM.VN() << "time:" << float(t1 - t0) / CLOCKS_PER_SEC << "s" << endl;        
        }

        //tri::io::Exporter<MyMesh>::Save(subM, pfileOut);
        vector<Point3f> result;
        for (auto vi = subM.vert.begin(); vi != subM.vert.end(); ++vi) if (!(*vi).IsD())
        {
            Point3f p_i;
            p_i[0] = (*vi).P()[0];
            p_i[1] = (*vi).P()[1];
            p_i[2] = (*vi).P()[2]; 
            result.push_back(p_i);
        }
        return result;    
    }

private:

    void Load_Resampling_Store() {
       
        std::vector<int> VertexId(m.vert.size());
        int index_R = 0;
        for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi) if (!(*vi).IsD())
        {
            VertexId[index_R] = index_R;
            index_R++;
            Point3f p_i;
            p_i[0] = (*vi).P()[0];
            p_i[1] = (*vi).P()[1];
            p_i[2] = (*vi).P()[2];
            pointOriginal.push_back(p_i);
            
        } 

        for (auto fi = m.face.begin(); fi != m.face.end(); ++fi) if (!(*fi).IsD()) {

            vector<int> face_i(3);
            face_i[0] = VertexId[tri::Index(m, (*fi).V(0))] ;//index of vertex per face   
            face_i[1] = VertexId[tri::Index(m, (*fi).V(1))] ;//index of vertex per face  
            face_i[2] = VertexId[tri::Index(m, (*fi).V(2))] ;//index of vertex per face  
            faceOriginal.push_back(face_i);

        }       
    
    }


};