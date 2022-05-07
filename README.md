# Adaptively-Isotropic-Remeshing

The code of adaptively isotropic remeshing.

If you want to use the original code, pls add public library below:

1. VCG libiary http://vcg.isti.cnr.it/vcglib/install.html
2. PCL 1.8+ https://pointclouds.org/downloads/
3. CGAL libiary 
4. Eigen3+
5. boost_1_67+

if you just want to output some results, you can use the exe file in EXE folder.

input parameters:
1. input file path (obj, ply)
2. if you want to subremesh and check face error in mesh (CAD mesh require),0: no; 1: yes
3. isotropic select; 1: isotropic; 2: adaptively isotropic 
4. accurate: 1 by default, more less, more accurate.

instance:

Data/Mesh/hand.obj 0 2 1
