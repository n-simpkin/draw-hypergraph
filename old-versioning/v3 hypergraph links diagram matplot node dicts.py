import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np
from scipy.interpolate import CubicSpline
from math import sqrt
from math import atan
from math import degrees
from math import cos
from math import sin
from math import radians

def setUpMatplotCanvas():
      fig,ax = plt.subplots()
      plt.xlim(-300,300)
      plt.ylim(-300,300)
      ax.set_aspect(1)
      return fig, ax
      
def drawNode(nodeCoords, ax, nodeRadius):
      nodeCircle = plt.Circle((nodeCoords[0],nodeCoords[1]), nodeRadius, fc="white", ec="black")
      ax.add_patch(nodeCircle)
      return ax
      
      
def setUpNodeDicts(nodeCoords):
      nodesInfo = []
      properties = ["coords","unitVector","adjustedUnitVector","associatedPolygonPoint"]
      keys = [None for i in range(5)]
      for i,node in enumerate(nodeCoords):
            nodesInfo.append(dict(zip(properties,keys)))
            nodesInfo[i]["coords"] = node
      return nodesInfo

def findCentroid(nodesInfo):
      centroid = [0,0]
      for node in nodesInfo:
            centroid[0] += node["coords"][0]
            centroid[1] += node["coords"][1]
      centroid[0] /= 3
      centroid[1] /= 3
      return centroid

def drawCentreCircle(centroid,radius,ax):
      #It looks a bit off but that's just how the circumcentre is.
      centreCircle = plt.Circle((centroid[0],centroid[1]),radius,fill=False)
      ax.add_patch(centreCircle)
      return ax

#find points of centre polygon, aligned correctly.
def constructCentrePolygon(nodesInfo,centroid,radius):
      polygonPoints = []
      for node in nodesInfo:
            vector = [node["coords"][i] - centroid[i] for i in range(2)] #minus node position from centroid for x and y. There are a few generators that look like this that will come up, they just do it for both the x and y.
            #print("vector of",node["coords"],"is",vector)
            vectorMag = sqrt((vector[0]**2) + (vector[1])**2) # Look out for negatives here
            node["unitVector"] = [vector[i]/vectorMag for i in range(2)]
            node["adjustedUnitVector"] = [node["unitVector"][i] + centroid[i] for i in range(2)] 
            #print("unit vector of",node["coords"],"is", node["unitVector"])
            polygonPoint = [(5*radius)*coord for coord in node["unitVector"]] #multiply unit vector by a multiple of the radius so the shape will sit outside the circle
            polygonPoint = [polygonPoint[i] + centroid[i] for i in range(2)] #Adjust the point so it's in relation to the centroid
            node["associatedPolygonPoint"] = polygonPoint
            #print("closestPolygonPoint of",node["coords"],"is", node["associatedPolygonPoint"])
            polygonPoints.append(polygonPoint)
      return polygonPoints

def constructCentrePolygonPointsProportion(nodesInfo,centroid,radius): #find points of centre polygon, aligned correctly.
      polygonPoints = []
      for node in nodesInfo:
            vector = [node["coords"][i] - centroid[i] for i in range(2)] #minus node position from centroid for x and y. 
            vectorMag = sqrt((vector[0]**2) + (vector[1])**2) 
            node["unitVector"] = [vector[i]/vectorMag for i in range(2)]
            polygonPoint = [(0.25)*coord for coord in vector] #multiply unit vector by a multiple of the radius so the shape will sit outside the circle
            polygonPoint = [polygonPoint[i] + centroid[i] for i in range(2)] #translate the point so its centred on the centroid rather than the origin
            node["associatedPolygonPoint"] = polygonPoint
            polygonPoints.append(polygonPoint)
      return polygonPoints

##def calcAngleOfRotation(polygonPoint,nextNodePolygonPoint):
##      #making a triangle of the centre polygon edge in question and the x axis, to find how far we need to rotate all the points so the two edge ones lie flat on the x axis.
##      trigWidth = polygonPoint[0] - nextNodePolygonPoint[0]
##      trigHeight = polygonPoint[1] - nextNodePolygonPoint[1]
##      #print("trig width",trigWidth,"trig height",trigHeight)
##      angle = 0
##      if trigHeight != 0:
##            angle = 90 + degrees(atan(trigWidth/trigHeight))
##            #print("angle",angle, "for", polygonPoint, "and", nextNodePolygonPoint)
##      return angle
##
##def rotatePoints(points, rotationPoint, rotationAngle):
##      #plt.plot(0,0,"bo")
##      rotationAngle = radians(rotationAngle)
##      rotatedPoints = []
##      for point in points:
##            translatedPoint = [point[i]-rotationPoint[i] for i in range(2)] #Translate points so the rotation point is at the origin.
##            rotatedX = (translatedPoint[0]*cos(rotationAngle)) - (translatedPoint[1]*sin(rotationAngle))
##            rotatedY = (translatedPoint[1]*cos(rotationAngle)) + (translatedPoint[0]*sin(rotationAngle))
##            rotatedX += rotationPoint[0] #Translate points back to where they were
##            rotatedY += rotationPoint[1]
##            #plt.plot(rotatedX,rotatedY,"yo")
##            rotatedPoints.append([rotatedX, rotatedY])
##      #print("rp", rotatedPoints)
##      return rotatedPoints
      
            
def drawCentrePolygonPoints(polygonPoints,plt):
      for point in polygonPoints:
            plt.plot(point[0], point[1], "ro")

def findCircleConnectionPoint(node1, node2, radius, centroid):
      direction = [node1["unitVector"][i] + node2["unitVector"][i] for i in range(2)]
      circlePoint = [radius*coord for coord in direction]
      circlePoint = [circlePoint[i] + centroid[i] for i in range(2)]
      plt.plot(circlePoint[0], circlePoint[1], "bo")
      return circlePoint
      
#draw the lines
def drawLines(nodesInfo, radius, centroid, ax):
      for currentNodeIndex in range(len(nodesInfo)):
      #for currentNodeIndex in range(1,2):
            nextNode = (currentNodeIndex + 1) % len(nodesInfo)
            circleConnectionPoint = findCircleConnectionPoint(nodesInfo[currentNodeIndex],nodesInfo[nextNode], radius, centroid)
            #rotationAngle = calcAngleOfRotation(nodesInfo[currentNodeIndex]["associatedPolygonPoint"],nodesInfo[nextNode]["associatedPolygonPoint"])
            pointsForParabola = [nodesInfo[currentNodeIndex]["associatedPolygonPoint"],circleConnectionPoint, nodesInfo[nextNode]["associatedPolygonPoint"]]
            #contrPoint = 10*(nodesInfo[currentNodeIndex]["adjustedUnitVector"] + nodesInfo[nextNode]["adjustedUnitVector"]) #can't just add silly
            #print([(nodesInfo[currentNodeIndex]["adjustedUnitVector"][i] + nodesInfo[nextNode]["adjustedUnitVector"][i]) for i in range(2)])
            contrPoint = [-(radius/2)*coord for coord in [(nodesInfo[currentNodeIndex]["unitVector"][i] + nodesInfo[nextNode]["unitVector"][i]) for i in range(2)]]
            contrPoint = [contrPoint[i] + centroid[i] for i in range(2)]
            print("cp",contrPoint)
            plt.plot(contrPoint[0],contrPoint[1],"yo")
            verts = [nodesInfo[currentNodeIndex]["associatedPolygonPoint"],contrPoint,nodesInfo[nextNode]["associatedPolygonPoint"]] #((-8.71- sqrt(17.5)),(68.33- sqrt(17.5))) ((-8.71+10),(78.24- 10))
            #middle val 10*-(sum of unit vectors)
            print(verts)
##            print("cur",nodesInfo[currentNodeIndex]["unitVector"])
##            print("nxt",nodesInfo[nextNode]["unitVector"],)
            #plt.axline(nodesInfo[currentNodeIndex]["coords"],nodesInfo[currentNodeIndex]["adjustedUnitVector"])
            #furthest circle point from node2 - 10? I think so. 10 is radius
            #(sum of unit vectors )back*10
            codes = [Path.MOVETO,Path.CURVE3, Path.CURVE3]
            bezier1 = patches.PathPatch(Path(verts, codes), fc="none")
            ax.add_patch(bezier1)
            plt.plot((nodesInfo[currentNodeIndex]["coords"][0],nodesInfo[currentNodeIndex]["associatedPolygonPoint"][0]), (nodesInfo[currentNodeIndex]["coords"][1],nodesInfo[currentNodeIndex]["associatedPolygonPoint"][1]), color="black")
##            if rotationAngle != 0:
##                  pointsForParabola = rotatePoints(pointsForParabola, nodesInfo[currentNodeIndex]["associatedPolygonPoint"], rotationAngle)
            #print("pp",pointsForParabola)
##            x = [point[0] for point in pointsForParabola]
##            y = [point[1] for point in pointsForParabola]
##            pairs = dict(zip(x,y))
##            x = sorted(x)
##            sortedY = []
##            for xval in x:
##                  sortedY.append(pairs[xval])
####            x = [nodesInfo[currentNodeIndex]["coords"][0],nodesInfo[currentNodeIndex]["associatedPolygonPoint"][0],circleConnectionPoint[0],nodesInfo[nextNode]["associatedPolygonPoint"][0], nodesInfo[nextNode]["coords"][0]]
####            y = [nodesInfo[currentNodeIndex]["coords"][1],nodesInfo[currentNodeIndex]["associatedPolygonPoint"][0], circleConnectionPoint[1], nodesInfo[nextNode]["associatedPolygonPoint"][0], nodesInfo[nextNode]["coords"][1]]
##            #print(x)
##            xvals = np.linspace(x,10)
##            spl = CubicSpline(x, sortedY, axis=13) # First generate spline function
##            #WONT ALWAYS BE A QUADRATIC AS NOT ALWAYS SYMMETRIC
##            y_smooth = spl(xvals) # then evalute for your interpolated points
##            #print([item.split() for item in xvals])
####            xvals = rotatePoints(xvals,pointsForParabola[0],rotationAngle)
####            y_smooth = rotatePoints(xvals,pointsForParabola[0],rotationAngle)
##            plt.plot(xvals, y_smooth)
                  
##            x = sorted([nodesInfo[currentNodeIndex]["coords"][0],circleConnectionPoint[0],nodesInfo[nextNode]["coords"][0]])
##            y = [nodesInfo[currentNodeIndex]["coords"][1],circleConnectionPoint[1], nodesInfo[nextNode]["coords"][1]]
            
##            x = [nodesInfo[currentNodeIndex]["associatedPolygonPoint"][0],circleConnectionPoint[0],nodesInfo[nextNode]["associatedPolygonPoint"][0]]
##            y = [nodesInfo[currentNodeIndex]["associatedPolygonPoint"][1],circleConnectionPoint[1], nodesInfo[nextNode]["associatedPolygonPoint"][1]]
##            pairs = dict(zip(x,y))
##            x = sorted(x)
##            sortedY = []
##            for xval in x:
##                  sortedY.append(pairs[xval])
####            x = [nodesInfo[currentNodeIndex]["coords"][0],nodesInfo[currentNodeIndex]["associatedPolygonPoint"][0],circleConnectionPoint[0],nodesInfo[nextNode]["associatedPolygonPoint"][0], nodesInfo[nextNode]["coords"][0]]
####            y = [nodesInfo[currentNodeIndex]["coords"][1],nodesInfo[currentNodeIndex]["associatedPolygonPoint"][0], circleConnectionPoint[1], nodesInfo[nextNode]["associatedPolygonPoint"][0], nodesInfo[nextNode]["coords"][1]]
##            xvals = np.linspace(x,10)
##            spl = CubicSpline(x, sortedY, axis=13) # First generate spline function
##            y_smooth = spl(xvals) # then evalute for your interpolated points
##            plt.plot(xvals, y_smooth)
##            #print("start of line located",nodesInfo[currentNodeIndex]["coords"], "end of line located",nodesInfo[currentNodeIndex]["associatedPolygonPoint"])
##            plt.plot((nodesInfo[currentNodeIndex]["coords"][0],nodesInfo[currentNodeIndex]["associatedPolygonPoint"][0]),(nodesInfo[currentNodeIndex]["coords"][1], nodesInfo[currentNodeIndex]["associatedPolygonPoint"][1]))
            #Focus, Directrix parabola




#Main program
nodesInfo = setUpNodeDicts([[-260,220],[260,220],[0,-220]]) #What if they're out of order? Nodes must be sequentially next to their spacewise neighbours.
#nodesInfo = setUpNodeDicts([[-260,-220],[260,-220],[0,220]])
radius = 25
nodeRadius = 30
#print(nodesInfo)

fig, ax = setUpMatplotCanvas()
for node in nodesInfo:
      ax = drawNode(node["coords"], ax, nodeRadius)
centroid = findCentroid(nodesInfo)
#centroid = [-50,-50]
drawCentreCircle(centroid,radius, ax)
polygonPoints = constructCentrePolygonPointsProportion(nodesInfo, centroid, radius)
drawCentrePolygonPoints(polygonPoints,plt)
drawLines(nodesInfo, radius, centroid, ax)
#drawLines([nodesInfo[1]],centroid,radius)
plt.show()
