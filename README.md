## Introduction

This project is the implementation of [CWF](https://github.com/Xrvitd/CWF) paper in Blender.

## Features

·Directly remesh the selected mesh object.

·Set target sampling value

## Todo
·Cleaning overlapping mesh and correcting surface normals

## Compile Manually

1. Install xmake

2. clone the project `git clone https://github.com/AIGODLIKE/Blender-CWF-Remesher`
   
3. cd to `Blender-CWF-Remesher/thirdparty/CWF`
   
4. Run config first
   
   1. for example macos with arm64 architecture `xmake f -p macosx -a arm64 -m release`

5. Run `xmake build -v remesh` to auto compile remesh module
   
   1. After build it will copy to Blender-CWF-Remesher/cxxlibs/youros/remesh_xxx_xxx.pyd
   
6. Install Blender-CWF-Remesher as blender addon
