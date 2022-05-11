   # -*- coding: utf-8 -*-
"""
Created on Tue May 11 09:11:24 2021

@author: DELL
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:34:14 2021

@author: DELL
"""



#   numpy.savetxt('xcenp'+repr(flag)+val+'.txt',xcenp,fmt='%f')
#   numpy.savetxt('ycenp'+repr(flag)+val+'.txt',ycenp,fmt='%f')
#   numpy.savetxt('elementnumber'+repr(flag)+val+'.txt', nel, fmt='%d')
#   numpy.savetxt('SCF'+repr(flag)+val+'.txt',scf,fmt='%f')
#   numpy.savetxt('nom_stress'+repr(flag)+val+'.txt',ns,fmt='%f')
    # numpy.savetxt('maxi_s11'+repr(flag)+val+'.txt',maxi_s11,fmt='%f')
    # numpy.savetxt('nodenumber6'+val+'.txt',non,fmt='%d')



# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.14-1 replay file
# Internal Version: 2014_06_05-03.41.02 134264s
# Run by user on Sat Oct 06 18:20:00 2018
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...


#length of the platelet = lth
#width of the platelet = wth
#aspect ratio of platelet =ar
#diameter of the hole =di
#position of the hole =psh
#thickness of matrix =thi
#thih=thi/2
#lthh=lth/2
#wthh=wth/2
#a11=X coordinates of center of the hole
#a12=Y coordinates of center of the hole
#w12 = width of the cutted portion
#di/w12 ratio = 0.1
#w11 = length of he cutted portion
import random
import numpy
# from collections import OrderedDict
from itertools import repeat
#import xlwt
import job
import os
from jobMessage import *


thi=5.0
di=10
lth=120
wth = 20.0
thih=thi/2
lthh=lth/2
wthh=wth/2
rad=di/2
u1 =[]
u2 = []
xx=[]
yy=[]
elst=[]
nel=[]
xcenp=[]
ycenp=[]
xcenpall=[]
ycenpall=[]
tot=[]
maxi_s11=[]
ns = []
scf = []
s11=[]
s22=[]
s33=[]
s44=[]

os.chdir(r"E:\\10D")
flag=1
n= 1
tim=0
k=0     
while flag<=n:
    con=1
    ok=0
    u1 =[]
    u2 = []
    xx=[]
    yy=[]
    e11 =[]
    e12 = []
    e22 = []
    from abaqus import *
    from abaqusConstants import *
    session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=121.633964538574, 
    height=118.363288879395)
    session.viewports['Viewport: 1'].makeCurrent()
    session.viewports['Viewport: 1'].maximize()
    from caeModules import *
    from driverUtils import executeOnCaeStartup
    executeOnCaeStartup()
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
    Mdb()


    lbls=[]
    session.viewports['Viewport: 1'].setValues(displayedObject=None)
    import sketch
    import part
    ###Whatever print functions we are using we do it to be shown in message area.
    print('iteration number '+repr(flag))

    a11=random.randint(300,330)      #last value wont count
    a12=random.randint(120,132)

    xcenpall.append(a11)
    ycenpall.append(a12)
    numpy.savetxt('a11.txt',xcenpall,fmt='%f')
    numpy.savetxt('a12.txt',ycenpall,fmt='%f')


    mdb.models.changeKey(fromName='Model-1', toName='Model-10')


    s = mdb.models['Model-10'].ConstrainedSketch(name='__profile__', 
    sheetSize=2000.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(thih, 0.0), point2=(lth+thih, wthh))
    s.rectangle(point1=(0.0, wthh+thi), point2=(lthh, wthh+thi+wth))
    s.rectangle(point1=(thih, wthh+thi+wth+thi), point2=(lth+thih, wthh+thi+wth+thi+wthh))
    s.rectangle(point1=(lthh+thi, wthh+thi), point2=(lth+thi, wthh+thi+wth))
    session.viewports['Viewport: 1'].view.fitView()
    p = mdb.models['Model-10'].Part(name='Part-1', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()

    #####PART2 ,8,9 for matrix

    s1 = mdb.models['Model-10'].ConstrainedSketch(name='__profile__', 
    sheetSize=2000.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)

    s1.Line(point1=(0.0, 0.0), point2=(thih, 0.0))
    s1.Line(point1=(thih, 0.0), point2=(thih, wthh))
    s1.Line(point1=(thih, wthh), point2=(lth+thih, wthh))
    s1.Line(point1=(lth+thih, wthh), point2=(lth+thih, 0.0))
    s1.Line(point1=(lth+thih, 0.0), point2=(lth+thi, 0.0))
    session.viewports['Viewport: 1'].view.fitView()
    s1.Line(point1=(lth+thi, 0.0), point2=(lth+thi, wthh+thi))
    s1.Line(point1=(lth+thi, wthh+thi), point2=(lthh+thi, wthh+thi))
    s1.Line(point1=(lthh+thi, wthh+thi), point2=(lthh+thi, wth+wthh+thi))
    session.viewports['Viewport: 1'].view.fitView()
    s1.Line(point1=(lthh+thi, wth+wthh+thi), point2=(lth+thi, wth+wthh+thi))
    s1.Line(point1=(lth+thi, wth+wthh+thi), point2=(lth+thi, wth+wth+thi+thi))
    session.viewports['Viewport: 1'].view.fitView()
    s1.Line(point1=(lth+thi, wth+wth+thi+thi), point2=(lth+thih, wthh+wthh+wth+thi+thi))
    s1.Line(point1=(lth+thih, wth+wth+thi+thi), point2=(lth+thih, wthh+wth+thi+thi))
    s1.Line(point1=(lth+thih, wthh+wth+thi+thi), point2=(thih, wthh+wth+thi+thi))
    s1.Line(point1=(thih, wthh+wth+thi+thi), point2=(thih, wthh+wthh+wth+thi+thi))
    s1.Line(point1=(thih, wthh+wthh+wth+thi+thi), point2=(0.0, wthh+wthh+wth+thi+thi))
    s1.Line(point1=(0.0, wthh+wthh+wth+thi+thi), point2=(0.0, wth+wthh+thi))
    s1.Line(point1=(0.0, wth+wthh+thi), point2=(lthh, wth+wthh+thi))
    s1.Line(point1=(lthh, wth+wthh+thi), point2=(lthh, wthh+thi))
    s1.Line(point1=(lthh, wthh+thi), point2=(0.0, wthh+thi))
    s1.Line(point1=(0.0, wthh+thi), point2=(0.0, 0.0))
    p = mdb.models['Model-10'].Part(name='Part-2', dimensionality=TWO_D_PLANAR, 
        type=DEFORMABLE_BODY)
    p.BaseShell(sketch=s1)
    s1.unsetPrimaryObject()

    mdb.models['Model-10'].Material(name='Material-1')
    mdb.models['Model-10'].materials['Material-1'].Elastic(table=((100.0, 0.3), ))
    mdb.models['Model-10'].Material(name='Material-2')
    mdb.models['Model-10'].materials['Material-2'].Elastic(table=((1000.0, 0.2), ))
    mdb.models['Model-10'].HomogeneousSolidSection(name='Section-1', 
        material='Material-1', thickness=None)
    mdb.models['Model-10'].HomogeneousSolidSection(name='Section-2', 
        material='Material-2', thickness=None)


    a = mdb.models['Model-10'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    a = mdb.models['Model-10'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-10'].parts['Part-1']
    a.Instance(name='Part-1-1', part=p, dependent=OFF)
    p = mdb.models['Model-10'].parts['Part-2']
    a.Instance(name='Part-2-1', part=p, dependent=OFF)
    a = mdb.models['Model-10'].rootAssembly
    a.InstanceFromBooleanMerge(name='Part-4', instances=(a.instances['Part-1-1'], 
        a.instances['Part-2-1'], ), keepIntersections=ON, 
        originalInstances=SUPPRESS, domain=GEOMETRY)
    a = mdb.models['Model-10'].rootAssembly
    a.makeIndependent(instances=(a.instances['Part-4-1'], ))
    a = mdb.models['Model-10'].rootAssembly

    p = mdb.models['Model-10'].parts['Part-4']
    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
        engineeringFeatures=ON)
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)

    platelet_face_point_1 = ((thih+lthh),wthh/2,0.0)
    platelet_face_point_2 = ((lthh/2),(wthh+thi+wthh),0.0)
    platelet_face_point_3 = ((lthh+thi+lthh/2),(wthh+thi+wthh),0.0)
    platelet_face_point_4 = ((thih+lthh),wth+wthh+thi+thi+wthh/2,0.0)
    platelet_face_1 = p.faces.findAt((platelet_face_point_1,))
    platelet_face_2 = p.faces.findAt((platelet_face_point_2,))
    platelet_face_3 = p.faces.findAt((platelet_face_point_3,))
    platelet_face_4 = p.faces.findAt((platelet_face_point_4,))
    platelet_region=(platelet_face_1, platelet_face_2, platelet_face_3, platelet_face_4,)
    p = mdb.models['Model-10'].parts['Part-4']
    p.SectionAssignment(region=platelet_region, sectionName='Section-2', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)


    p = mdb.models['Model-10'].parts['Part-4']
    f = p.faces
    matrix_face_point_1 = ((thih/2),(wthh/2),0.0)
    matrix_face_point_2 = ((lth+thih+(thih/2)),(wthh/2),0.0)
    matrix_face_point_3 = (lthh+thih,wthh+thi+wthh,0.0)
    matrix_face_point_4 = ((lthh+thih)/2,wthh+(thi/2),0.0)
    matrix_face_point_5 = (lthh+thih+(lthh+thih)/2,wthh+(thi/2),0.0)
    matrix_face_point_6 = ((lthh+thih)/2,wth+thi+(thi/2)+wthh,0.0)
    matrix_face_point_7 = (lthh+thih+(lthh+thih)/2,wth+thi+(thi/2)+wthh,0.0)
    matrix_face_point_8 = ((thih/2),(wthh/2)+thi+thi+wthh+wth,0.0)
    matrix_face_point_9 = ((lth+thih+(thih/2)),(wthh/2)+thi+thi+wthh+wth,0.0)
    matrix_face_1 = p.faces.findAt((matrix_face_point_1,))
    matrix_face_2 = p.faces.findAt((matrix_face_point_2,))
    matrix_face_3 = p.faces.findAt((matrix_face_point_3,))
    matrix_face_4 = p.faces.findAt((matrix_face_point_4,))
    matrix_face_5 = p.faces.findAt((matrix_face_point_5,))
    matrix_face_6 = p.faces.findAt((matrix_face_point_6,))
    matrix_face_7 = p.faces.findAt((matrix_face_point_7,))
    matrix_face_8 = p.faces.findAt((matrix_face_point_8,))
    matrix_face_9 = p.faces.findAt((matrix_face_point_9,))
    matrix_region=(matrix_face_1,matrix_face_2,matrix_face_3,matrix_face_4,matrix_face_5,matrix_face_6,matrix_face_7,
        matrix_face_8,matrix_face_9,)
    p = mdb.models['Model-10'].parts['Part-4']
    p.SectionAssignment(region=matrix_region, sectionName='Section-1', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)



    s = mdb.models['Model-10'].ConstrainedSketch(name='__profile__', 
    sheetSize=2000.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s1=a11-150
    s2=a11+150
    s3=a12-60
    s4=a12+60
    s.rectangle(point1=(s1,s3), point2=(s2,s4))
    s.CircleByCenterPerimeter(center=(a11,a12), point1=(a11+rad,a12))
    p = mdb.models['Model-10'].Part(name='Part-3', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
    p = mdb.models['Model-10'].parts['Part-3']
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()

    for i in range(0,100):
        p1 = mdb.models['Model-10'].parts['Part-4']
        session.viewports['Viewport: 1'].setValues(displayedObject=p1)
        p = mdb.models['Model-10'].Part(name='rve'+repr(i), 
        objectToCopy=mdb.models['Model-10'].parts['Part-4'])
        session.viewports['Viewport: 1'].setValues(displayedObject=p)
        p = mdb.models['Model-10'].parts['rve'+repr(i)]
        a = mdb.models['Model-10'].rootAssembly
        a.Instance(name='rve'+repr(i)+'-1', part=p, dependent=OFF)
    a = mdb.models['Model-10'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    a.features['Part-4-1'].suppress()
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    dis1=lth+thih+thih
    dis2=0.0
    b1=0
    for j in range(1,11):
      for i in range(0,10):
        a1 = mdb.models['Model-10'].rootAssembly
        a1.translate(instanceList=('rve'+repr(b1)+'-1', ), vector=(i*dis1, dis2, 0.0))
        b1=b1+1
      dis2=j*(wth+wth+thi+thi)   
    session.viewports['Viewport: 1'].view.fitView()
    a = mdb.models['Model-10'].rootAssembly
    q=[0]*len(range(0,100))
    for i in range(0,100):
     q[i]=a.instances['rve'+repr(i)+'-1']
    a.InstanceFromBooleanMerge(name='Part-5', instances=q, 
    originalInstances=SUPPRESS, domain=GEOMETRY)
    a.makeIndependent(instances=(a.instances['Part-5-1'], ))
    p = mdb.models['Model-10'].parts['Part-5']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)


    p = mdb.models['Model-10'].parts['Part-3']
    a.Instance(name='Part-3-1', part=p, dependent=OFF)
    p = mdb.models['Model-10'].parts['Part-5']
    a.Instance(name='Part-5-1', part=p, dependent=OFF)
    a = mdb.models['Model-10'].rootAssembly
    a.InstanceFromBooleanCut(name='Part-6',
    instanceToBeCut=mdb.models['Model-10'].rootAssembly.instances['Part-5-1'],
    cuttingInstances=(a.instances['Part-3-1'], ), originalInstances=DELETE)
    p1 = mdb.models['Model-10'].parts['Part-6']
    a.makeIndependent(instances=(a.instances['Part-6-1'], ))
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)

    a = mdb.models['Model-10'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    a1 = mdb.models['Model-10'].rootAssembly
    p = mdb.models['Model-10'].parts['Part-5']
    a1.Instance(name='Part-5-1', part=p, dependent=OFF)
    p = mdb.models['Model-10'].parts['Part-6']
    a1.Instance(name='Part-6-1', part=p, dependent=OFF)
    a1 = mdb.models['Model-10'].rootAssembly
    a1.InstanceFromBooleanCut(name='Part-7', 
    instanceToBeCut=mdb.models['Model-10'].rootAssembly.instances['Part-5-1'], 
    cuttingInstances=(a1.instances['Part-6-1'], ), originalInstances=DELETE)
    a.makeIndependent(instances=(a.instances['Part-7-1'], ))
    p1 = mdb.models['Model-10'].parts['Part-7']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)

    s = mdb.models['Model-10'].ConstrainedSketch(name='__profile__', 
    sheetSize=2000.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s8=a11
    s9=a12
    s6=a11+ 12.5*di
    s7=a12+5*di
    s.rectangle(point1=(s8,s9), point2=(s6,s7))
    #s.ArcByCenterEnds(center=(a11, a12), point1=(a11, a12 + rad), point2=(a11+rad, a12), direction=CLOCKWISE)
    p = mdb.models['Model-10'].Part(name='Part-9', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
    p = mdb.models['Model-10'].parts['Part-9']
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()

    s = mdb.models['Model-10'].ConstrainedSketch(name='__profile__', 
    sheetSize=2000.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.ArcByCenterEnds(center=(a11, a12), point1=(a11, a12 + rad), point2=(a11+rad, a12), direction=CLOCKWISE)
    s.Line(point1=(a11, a12), point2=(a11, a12 + rad))
    s.Line(point1=(a11, a12), point2= (a11+rad, a12))
    p = mdb.models['Model-10'].Part(name='Part-8', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
    p = mdb.models['Model-10'].parts['Part-8']
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()

    p = mdb.models['Model-10'].parts['Part-8']
    a.Instance(name='Part-8-1', part=p, dependent=OFF)
    p = mdb.models['Model-10'].parts['Part-9']
    a.Instance(name='Part-9-1', part=p, dependent=OFF)
    a = mdb.models['Model-10'].rootAssembly
    a.InstanceFromBooleanCut(name='Part-10',
    instanceToBeCut=mdb.models['Model-10'].rootAssembly.instances['Part-9-1'],
    cuttingInstances=(a.instances['Part-8-1'], ), originalInstances=DELETE)
    p1 = mdb.models['Model-10'].parts['Part-10']
    a.makeIndependent(instances=(a.instances['Part-10-1'], ))
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)

    p = mdb.models['Model-10'].parts['Part-10']
    a.Instance(name='Part-10-1', part=p, dependent=OFF)
    p = mdb.models['Model-10'].parts['Part-7']
    a.Instance(name='Part-7-1', part=p, dependent=OFF)
    a = mdb.models['Model-10'].rootAssembly
    a.InstanceFromBooleanCut(name='Part-11',
    instanceToBeCut=mdb.models['Model-10'].rootAssembly.instances['Part-7-1'],
    cuttingInstances=(a.instances['Part-10-1'], ), originalInstances=DELETE)
    p1 = mdb.models['Model-10'].parts['Part-11']
    a.makeIndependent(instances=(a.instances['Part-11-1'], ))
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)

    a = mdb.models['Model-10'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    a1 = mdb.models['Model-10'].rootAssembly
    p = mdb.models['Model-10'].parts['Part-7']
    a1.Instance(name='Part-7-1', part=p, dependent=OFF)
    p = mdb.models['Model-10'].parts['Part-11']
    a1.Instance(name='Part-11-1', part=p, dependent=OFF)
    a1 = mdb.models['Model-10'].rootAssembly
    a1.InstanceFromBooleanCut(name='Part-12', 
    instanceToBeCut=mdb.models['Model-10'].rootAssembly.instances['Part-7-1'], 
    cuttingInstances=(a1.instances['Part-11-1'], ), originalInstances=DELETE)
    a.makeIndependent(instances=(a.instances['Part-12-1'], ))
    p1 = mdb.models['Model-10'].parts['Part-12']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    ##cutting of the model from big


    #mdb.models['Model-10'].StaticStep(name='Step-1', previous='Initial')
    #session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')

    mdb.models['Model-10'].StaticStep(name='Step-1', previous='Initial', 
    timePeriod=0.1, initialInc=0.1, minInc=1e-06, maxInc=0.1)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)

    mdb.models['Model-10'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U','UT','UR', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'COORD'))


    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
    a = mdb.models['Model-10'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    f1 = a.instances['Part-12-1'].faces
    t = a.MakeSketchTransform(sketchPlane=f1[0], sketchPlaneSide=SIDE1, origin=(
    0, 0, 0.0))
    s = mdb.models['Model-10'].ConstrainedSketch(name='__profile__', 
    sheetSize=2000.0, gridSpacing=5.00, transform=t)
    s.setPrimaryObject(option=SUPERIMPOSE)
    a = mdb.models['Model-10'].rootAssembly
    a.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.rectangle(point1=(a11-4*rad, a12-4*rad), point2=(a11+4*rad, a12+4*rad))
    #s.rectangle(point1=(a11-1.8*rad, a12-1.8*rad), point2=(a11+1.8*rad, a12+1.8*rad))
    s.rectangle(point1=(s8+.5, s9+.5), point2=(s6-.5,s7-.5 ))
    #s.rectangle(point1=(a11-3*rad, a12-3*rad), point2=(a11+3*rad, a12+3*rad))
    s.rectangle(point1=(a11-2*rad, a12-2*rad), point2=(a11+2*rad, a12+2*rad))
    #s.Line(point1=(a11, a12-20), point2=(a11, a12 + 20))
    a.PartitionFaceBySketch(faces=f1, sketch=s)
    s.unsetPrimaryObject()



    a = mdb.models['Model-10'].rootAssembly
    e1 = a.instances['Part-12-1'].edges
    bon = e1.getByBoundingBox(s8, s9, 0.0, s8, s7, 0.0)+e1.getByBoundingBox(s6, s9, 0.0, s6, s7, 0.0) \
    +e1.getByBoundingBox(s8, s7, 0.0, s6, s7, 0.0)+e1.getByBoundingBox(s8, s9, 0.0, s6, s9, 0.0) \
    +e1.getByBoundingBox(s8+.5, s9+.5, 0.0, s8+.5, s7-.5, 0.0)+e1.getByBoundingBox(s6-.5, s9+.5, 0.0, s6-.5, s7-.5, 0.0)\
    +e1.getByBoundingBox(s8+.5, s7-.5, 0.0, s6-.5, s7-.5, 0.0)+e1.getByBoundingBox(s8+.5, s9+.5, 0.0, s6-.5, s9+.5, 0.0)

    a.Set(edges=bon, name='boundaries')


    a = mdb.models['Model-10'].rootAssembly
    f1 = a.instances['Part-12-1'].faces
    e1 = a.instances['Part-12-1'].edges
    Face=f1.getByBoundingBox(s8, s9, 0.0, s6, s7, 0.0)
    a.Set(faces=Face, name='cut6')
    region=a.sets['cut6']
    a = mdb.models['Model-10'].rootAssembly.sets['cut6'].faces  
    b = mdb.models['Model-10'].rootAssembly.sets['cut6']
    elemType1 = mesh.ElemType(elemCode=CPE4,elemLibrary=STANDARD)
    a = mdb.models['Model-10'].rootAssembly
    a.setElementType(regions=b, elemTypes=(elemType1,))

    #face1 = f1.getByBoundingBox(a11-4*rad, a12-4*rad, 0.0, a11+4*rad, a12+4*rad, 0.0)
    #a.Set(faces=face1, name='reqreg')
    a = mdb.models['Model-10'].rootAssembly
    f1 = a.instances['Part-12-1'].faces
    e1 = a.instances['Part-12-1'].edges
    Face=f1.getByBoundingBox(a11-4*rad, a12-4*rad, 0.0, a11+4*rad, a12+4*rad, 0.0)
    a.Set(faces=Face, name='cut66')
    region=a.sets['cut66']
    a = mdb.models['Model-10'].rootAssembly.sets['cut66'].faces  
    b = mdb.models['Model-10'].rootAssembly.sets['cut66']
    elemType1 = mesh.ElemType(elemCode=CPE4,elemLibrary=STANDARD)
    a = mdb.models['Model-10'].rootAssembly
    a.setElementType(regions=b, elemTypes=(elemType1,))

    remesh3 = e1.getByBoundingSphere((a11, a12, 0.0) , 1.01*rad)
    a.Set(edges=remesh3, name='circle')
    a = mdb.models['Model-10'].rootAssembly
    e1 = a.instances['Part-12-1'].edges
    remesh2 = e1.getByBoundingBox(a11-2*rad, a12-2*rad, 0.0, a11+2*rad, a12+2*rad, 0.0)
    a.Set(edges=remesh2, name='smallbox')
    a = mdb.models['Model-10'].rootAssembly
    e1 = a.instances['Part-12-1'].edges
    remesh1 = e1.getByBoundingBox(a11-4*rad, a12-4*rad, 0.0, a11+4*rad, a12+4*rad, 0.0)
    a.Set(edges=remesh1, name='bigbox')
    #remesh4 = e1.getByBoundingBox(a11-1.8*rad, a12-1.8*rad, 0.0, a11+1.8*rad, a12 + 1.8*rad, 0.0)
    #a.Set(edges=remesh2, name='smallbox')
    #a = mdb.models['Model-10'].rootAssembly
    #e1 = a.instances['Part-7-1'].edges

    a.seedEdgeBySize(edges=bon, size=1, deviationFactor=0.1, 
    constraint=FINER)
    a.seedEdgeBySize(edges=remesh1, size=1, deviationFactor=0.1, 
    constraint=FINER)
    a.seedEdgeBySize(edges=remesh2, size= 0.15, deviationFactor=0.1, 
    constraint=FINER)
    a.seedEdgeBySize(edges=remesh3, size= 0.15, deviationFactor=0.1, 
    constraint=FINER)
    #a.seedEdgeBySize(edges=remesh4, size= 0.035, deviationFactor=0.1, 
    # constraint=FINER)



    a = mdb.models['Model-10'].rootAssembly  
    partInstances =(a.instances['Part-12-1'], ) 
    a.seedPartInstance(regions=partInstances, size=3, deviationFactor=0.1, 
    minSizeFactor=0.1)



    a.generateMesh(regions=partInstances)

    a = mdb.models['Model-10'].rootAssembly
    e1 = a.instances['Part-12-1'].edges



    be = e1.getByBoundingBox(s8, s9, 0.0, s6, s9, 0.0)
    te = e1.getByBoundingBox(s8, s7, 0.0, s6, s7, 0.0)
    le = e1.getByBoundingBox(s8, s9, 0.0, s8, s7, 0.0)
    re = e1.getByBoundingBox(s6, s9, 0.0, s6, s7, 0.0)
    a.Set(edges=be, name='be')
    a.Set(edges=te, name='te')
    a.Set(edges=le, name='le')
    a.Set(edges=re, name='re')
    a2 = mdb.models['Model-10'].rootAssembly

    # region = a2.sets['te']
    # mdb.models['Model-10'].DisplacementBC(name='TE', createStepName='Step-1', 
    # region=region, u1=-0.1*50*flag, u2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
    # distributionType=UNIFORM, fieldName='', localCsys=None)
    # region = a2.sets['le']
    # mdb.models['Model-10'].DisplacementBC(name='LE', createStepName='Step-1', 
    # region=region, u1=-0.015, u2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
    # distributionType=UNIFORM, fieldName='', localCsys=None)
    region = a2.sets['re']
    mdb.models['Model-10'].DisplacementBC(name='RE', createStepName='Step-1', 
    region=region, u1=0.015, u2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
    distributionType=UNIFORM, fieldName='', localCsys=None)

    region = a2.sets['be'] 
    mdb.models['Model-10'].YsymmBC(name='BC-1', createStepName='Initial', region=region, localCsys=None)
    # a = mdb.models['Model-1'].rootAssembly
    # e1 = a.instances['Part-1-1'].edges

    region = a2.sets['le']
    mdb.models['Model-10'].XsymmBC(name='BC-2', createStepName='Initial', region=region, localCsys=None)


    # region = a.sets['re']
    # mdb.models['Model-10'].XsymmBC(name='RE', createStepName='Step-1', 
    # region=region, localCsys=None)


    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, loads=OFF, 
    bcs=OFF, predefinedFields=OFF, connectors=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)

    #   # mdb.Job(name='Jobinput', model='Model-10', description='', type=ANALYSIS, 
    #   #     atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    #   #     memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    #   #     explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    #   #     modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    #   #     scratch='', multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    #   # mdb.jobs['Jobinput'].writeInput(consistencyChecking=OFF)



    mdb.Job(name='simulatid', model='Model-10', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', multiprocessingMode=DEFAULT, numCpus=4,numDomains=4 , numGPUs=0)
    mdb.jobs['simulatid'].submit(consistencyChecking=OFF)

    mdb.jobs['simulatid'].waitForCompletion()


    job_name='simulatid'


    sta=mdb.jobs[job_name].status
    st=str(sta)



    odbPath='E:\\10D\\'+job_name+'.odb'                   ##Path
    odb=session.openOdb(odbPath)
    session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(deformationScaling=UNIFORM, uniformScaleFactor=0)
    lastFrame=odb.steps['Step-1'].frames[-1]
    disp=lastFrame.fieldOutputs['U']                      ##S-stress ##U-displacement
    reg=odb.rootAssembly.nodeSets['CUT66']                 ##Set name in capital ##only working for element sets  
    asd=disp.getSubset(region=reg,position=NODAL)
    locxo=s1
    locyo=s3
    for do in asd.values:   
        ff=do.nodeLabel
        u1r = do.data[0]
        u2r = do.data[1] 
        u1.append(u1r)
        u2.append(u2r)  



    pos=lastFrame.fieldOutputs['COORD']                      ##S-stress ##U-displacement
    reg=odb.rootAssembly.nodeSets['CUT66']                 ##Set name in capital ##only working for element and node sets  
    asp=pos.getSubset(region=reg,position=NODAL)
    for do in asp.values:   
        ff=do.nodeLabel
        xxr = do.data[0] 
        yyr = do.data[1] 
        locxo=s1
        locyo=s3
        xr=xxr-locxo
        yr=yyr-locyo 
        xx.append(xr)
        yy.append(yr)





    numpy.savetxt('x'+repr(flag)+'.txt',xx,fmt='%f')
    numpy.savetxt('y'+repr(flag)+'.txt',yy,fmt='%f')
    numpy.savetxt('u1'+repr(flag)+'.txt',u1,fmt='%f')
    numpy.savetxt('u2'+repr(flag)+'.txt',u2,fmt='%f')
    print(st)



    # #     # if st =='COMPLETED':
    # #     #   tim=0
    # #     #   ok=1
    # #     #   os.remove('Job.inp')

    # #     # else:
    # #     #   sys.exc_clear()
    # #     #   tim=1
    # #     #   os.remove(job_name+'.lck')
    # #     # #     k=k+1
    # #     #   os.remove('Job.inp')

    # #     # if ok==1:

    #odb = session.odbs['E:\10D\\Ar20\\'+job_name+'.odb']


    # #                   ####To store number of elements
    # #     ##To get fieldoutputs(stress,diplacement etc..) without writing xydata


    # # session.Path(name='Path-2', type=POINT_LIST, expression=((a11 , a12-4*rad, 0.0), (
    # #     a11 , a12-rad, 0.0)))
    # # session.Path(name='Path-1', type=POINT_LIST, expression=((a11, a12+rad, 0.0), (
    # #     a11 , a12+4*rad, 0.0)))
    # # # session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    # # #   variableLabel='E', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    # # #   'E22'))
    # # # pth = session.paths['Path-1']



    # # # session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=True, 
    # # #    pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
    # # # #     x0 = session.xyDataObjects['XYData-1']


    # # #   total=0
    # # #   for i in range(0,(len(x0)-1)):
    # # #       a1=list(x0[i])
    # # #       a2=list(x0[i+1])
    # # #       c=.5*(a2[0]-a1[0])*(a1[1]+a2[1])
    # # #       total=c+total
    # # #   nom_stress = total/(a2[0]-di)

    # # #   ns.append(nom_stress)
    # # #   scf.append(maxi_s11/nom_stress)

    # # #   save('AR20')
    flag+=1