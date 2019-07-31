# raytracing
Simple raytracing program with BVH and threading 

Here are some sample images produced by this program.

<img src=https://github.com/gkcheong/raytracing/blob/master/Images/arm.bmp width="400" height="400" /> |
<img src=https://github.com/gkcheong/raytracing/blob/master/Images/arm-top.bmp width="400" height="400" />
-------------------------------------------------------------------------------------------------------------------
<img src=https://github.com/gkcheong/raytracing/blob/master/Images/dragon_side_new.bmp width="400" height="400" /> |
<img src=https://github.com/gkcheong/raytracing/blob/master/Images/gear.bmp width="400" height="400" />


For some of the largest input files the usage of BVH and threading can give up to 1000x 
acceleration over simple compiler optimization (-O3)

The input file format is specific for this program and can be found in input/SceneFile.pdf
along with some input files for testing

This program is made for my auditing of CSCI 5607 in the University of Minnesota
