## Introduction

This project is the implementation of [CWF](https://github.com/Xrvitd/CWF) paper in Blender(3.6+).

## Features

·Directly remesh the selected mesh object.

·Set target sampling value

·Using pre simplification to improve computing speed

·Multi-Language Support(EN/zh-CN)

![image](https://github.com/AIGODLIKE/Blender-CWF-Remesher/assets/116185401/27dd2759-2ca2-4543-a909-cc5554bcdf37)

## Tips
1. Avoid [non-manifold ](https://docs.blender.org/manual/en/4.1/glossary/index.html#term-Non-manifold)

·Some borders and holes (edges with only a single connected face), as faces have no thickness.

·Edges and vertices not belonging to any face (wire).

·Edges connected to three or more faces (interior faces).

·Vertices belonging to faces that are not adjoining (e.g. two cones sharing the vertex at the apex).

2. Generate number of faces through sampling control（Usually, the larger the sampling, the more faces there are）

3. Simplification can improve speed（The smaller the simplification value, the faster it will be, but it will reduce the effectiveness of the CWF algorithm）


## Compile Manually

1. Install xmake

2. clone the project `git clone https://github.com/AIGODLIKE/Blender-CWF-Remesher`
   
3. cd to `Blender-CWF-Remesher/thirdparty/CWF`
   
4. Run config first
   
   1. for example macos with arm64 architecture `xmake f -p macosx -a arm64 -m release`

5. Run `xmake build -v remesh` to auto compile remesh module
   
   1. After build it will copy to Blender-CWF-Remesher/cxxlibs/youros/remesh_xxx_xxx.pyd
   
6. Install Blender-CWF-Remesher as blender addon
