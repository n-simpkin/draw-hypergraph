import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np #rewrite all wiht numpy arrays, eliminates all my clunky two generators, much more semantic.
from scipy.interpolate import CubicSpline
from math import sqrt
from math import atan
from math import degrees
from math import cos
from math import sin
from math import radians
from math import pi

'''
TODOS
Stop the overlapping - does this require a change to my whole approach?
The straight lines don't always flow in.
Hypergraph class for holding the information
How to avoid overlaps? Hold eges and thier circumcentres?
'''


def setUpMatplotCanvas():
      fig,ax = plt.subplots()
      plt.xlim(-300,300)
      plt.ylim(-300,300)
      ax.set_aspect(1)
      return fig, ax
      
def setUpNodeDicts(nodeCoords):
      nodesInfo = []
      properties = ["coords","unitVector","associatedPolygonPoint","circleConnectionPointDirection"]
      keys = [None for i in range(4)]
      for i,node in enumerate(nodeCoords):
            nodesInfo.append(dict(zip(properties,keys)))
            nodesInfo[i]["coords"] = node
      return nodesInfo

def findCentroid(nodesInfo):
      centroid = [0,0]
      #Take mean average of all coordinates to find the centroid
      for node in nodesInfo: 
            centroid = [centroid[i] + node["coords"][i] for i in range(2)] #There are a few generators that look like this that come up, they just perform the operation for both the x and y.
      centroid = [centroid[i]/3 for i in range(2)]
      return centroid

def calculateCentrePolygonPointsProportion(nodesInfo,centroid,radius, middleConnectorSize): #find points of centre polygon, aligned correctly.
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

def findCircleConnectionPointProportion(node1, node2, radius, centroid): #Should be fixed
      #The curved line must connect to the point on the circle that is equidistant from the two node points (the closest one as there will be two).
      #The direction of this point on the circle from the centroid can be found by adding the unit vectors (the direction of the node from the centroid) of the two nodes together.
      #Second above no longer always true as circle connecting point not in the middle of those directions.

      #Find point directly in the middle
      resultantVector = [node1["unitVector"][i] + node2["unitVector"][i] for i in range(2)] #Doesn't produce a unit vector so the point isn't actually radius distance away.

      
##      resultantVector = [node1["associatedPolygonPoint"][i] + node2["associatedPolygonPoint"][i] for i in range(2)]
##      resultantVector = [resultantVector[i] - centroid[i] for i in range(2)] #Put the vector in relation to 0,0 as everything else works in relation to 0,0
      
      #plt.plot(resultantVector[0],resultantVector[1],"bo")
      #plt.plot(centroid[0],centroid[1],"yo")
      #plt.plot(0,0,"yo")
      
      #Normalise
      vectorMag = sqrt((resultantVector[0]**2) + (resultantVector[1])**2) 
      direction = [resultantVector[i]/vectorMag for i in range(2)]
      #node1["circleConnectionPointDirection"] = direction
      
##      direction = [direction[i] + centroid[i] for i in range(2)]
      #plt.plot(direction[0],direction[1],"go")
##      direction = [direction[i] - centroid[i] for i in range(2)]
      
      circlePoint = [(radius*coord) for coord in direction] #The multiply the direction (of length one) by the circle radius to find the point on the circle needed.
##      a = [node1["unitVector"][i]*radius for i in range(2)]
##      b = [node2["unitVector"][i]*radius for i in range(2)]
##      circlePoint = [a[i] + b[i] for i in range(2)]
      circlePoint = [circlePoint[i] + centroid[i] for i in range(2)] #All this will have been performed with the centre as 0,0, so translate the point so its in relation to the centorid
      plt.plot(circlePoint[0],circlePoint[1],"yo")
      return circlePoint

def mag(coords):
            return sqrt((coords[0]**2) + (coords[1]**2))
      
def calcCentreCurveVertsProportion(nodeFrom,nodeTo, radius, centroid):
      #plt.plot((nodeFrom["associatedPolygonPoint"][0],nodeTo["associatedPolygonPoint"][0]),(nodeFrom["associatedPolygonPoint"][1], nodeFrom["associatedPolygonPoint"][1]))
      circleConnectionPoint = findCircleConnectionPointProportion(nodeFrom,nodeTo, radius, centroid)
      #plt.plot(circleConnectionPoint[0], circleConnectionPoint[1], "bo")
##      resultantVector = [(nodeFrom["unitVector"][i] + nodeTo["unitVector"][i]) for i in range(2)]
##      vectorMag = sqrt((resultantVector[0]**2) + (resultantVector[1])**2)
##      controlPointDirection = [resultantVector[i]/vectorMag for i in range(2)]
      #controlPointDirection = nodeFrom["circleConnectionPointDirection"]
      #controlPointDirection = [controlPointDirection[i] - centroid[i] for i in range(2)]
      #Control point definitley always has to be directly between the points. 
      #Turning point of curve is always half the distance from the start points to the control point.
      #Therefore control point must be 2* the distance from the polygon points to the circle point, in the direction of the added node Unit vectors.
##      vector = [circleConnectionPoint[i] - nodeFrom["associatedPolygonPoint"][i] for i in range(2)]
##      vector = [vector[i] + centroid[i] for i in range(2)]
##      vector = [vector[i] - centroid[i] for i in range(2)]
##      print(vector)
      #plt.plot(controlPointDirection[0],controlPointDirection[1],"bo")
      #midpointOfPolygonPoints = [(nodeFrom["associatedPolygonPoint"][i] + nodeTo["associatedPolygonPoint"][i])/2 for i in range(2)]
      #Midpoint no longer relevant as tp of curve not alwyas in the middle.

##      controlPointDirection = [nodeFrom["associatedPolygonPoint"][i] + nodeTo["associatedPolygonPoint"][i] for i in range(2)]
##      controlPointDirection = [controlPointDirection[i] / mag(controlPointDirection) for i in range(2)]
##      print(controlPointDirection)
      #controlPointDirection = [controlPointDirection[i] + centroid[i] for i in range(2)] #Put the vector in relation to 0,0 as everything else works in relation to 0,0
      #print(controlPointDirection)
      #plt.plot(controlPointDirection[0],controlPointDirection[1], "bo") 

      magNodeTo = mag(nodeTo["associatedPolygonPoint"])
      magNodeFrom = mag(nodeFrom["associatedPolygonPoint"])
      print("mnf", magNodeFrom)
      print("mnt", magNodeTo)
      
      if magNodeTo >= magNodeFrom: #What about if equal?
            print("1st")
            proportionOfPointOnLine = magNodeTo / magNodeFrom
            vector = [nodeFrom["associatedPolygonPoint"][i] - nodeTo["associatedPolygonPoint"][i] for i in range(2)]
            pointOnLine = [(vector[i] / (proportionOfPointOnLine+1)) for i in range(2)]
            pointOnLine = [pointOnLine[i] + nodeTo["associatedPolygonPoint"][i] for i in range(2)] # adjust to being from pointTo
      else:
            print("2nd")
            proportionOfPointOnLine = magNodeFrom / magNodeTo
            print("proportion",proportionOfPointOnLine)
            vector = [nodeTo["associatedPolygonPoint"][i] - nodeFrom["associatedPolygonPoint"][i] for i in range(2)]
            print(vector)
            pointOnLine = [(vector[i] / (proportionOfPointOnLine+1)) for i in range(2)]
            pointOnLine = [pointOnLine[i] + nodeFrom["associatedPolygonPoint"][i] for i in range(2)] # adjust to being from pointTo

        
      #linePoint =
      #pointOnLine = [-72,18]
      plt.plot(pointOnLine[0],pointOnLine[1], "go")

      controlPointDirection = [pointOnLine[i] - circleConnectionPoint[i] for i in range(2)]
      controlPointDirection = [controlPointDirection[i]/mag(controlPointDirection) for i in range(2)]

      controlPointDirection = [nodeFrom["unitVector"][i] + nodeTo["unitVector"][i] for i in range(2)]
      controlPointDirection = [controlPointDirection[i]/mag(controlPointDirection) for i in range(2)]
      
      pointOnLine = [(nodeFrom["associatedPolygonPoint"][i] + nodeTo["associatedPolygonPoint"][i])/2 for i in range(2)] #Make it midpoint
      vector = [pointOnLine[i] - circleConnectionPoint[i] for i in range(2)] #???
      vectorMag = sqrt((vector[0]**2) + (vector[1]**2))
      
      print("vm",vectorMag)
      print("ctrl pt direction", controlPointDirection)
      controlPoint = [vectorMag*(-coord)*2 for coord in controlPointDirection] # I don't think I can just use 2 anymore actually idk
      #print("ctrl pt", controlPoint)
      #controlPoint = [controlPoint[i] + pointOnLine[i] for i in range(2)]
      controlPoint = [controlPoint[i] + pointOnLine[i] for i in range(2)]
      #pFrom = (0,122.46)
      #controlPoint = [controlPoint[i] + centroid[i] for i in range(2)] #Here is the problem.
      #controlPoint = [controlPoint[i] + pFrom[i] for i in range(2)]
      #plt.plot(controlPoint[0],controlPoint[1], "yo")
      #controlPoint = [-90,-57]
      plt.plot(controlPoint[0],controlPoint[1], "bo") 
      curveVerts = [nodeFrom["associatedPolygonPoint"],controlPoint,nodeTo["associatedPolygonPoint"]]
      return curveVerts

def drawCentreCurve(curveVerts):
      codes = [Path.MOVETO,Path.CURVE3, Path.CURVE3]
      bezier1 = patches.PathPatch(Path(curveVerts, codes), fc="none")
      ax.add_patch(bezier1)

def drawNodeToPolygonLine(nodeCoords, polygonPointCoords):
      plt.plot((nodeCoords[0], polygonPointCoords[0]), (nodeCoords[1], polygonPointCoords[1]), color="black")

def drawCentreCircle(centroid,radius,ax):
      #It looks a bit off but that's just how the circumcentre is.
      centreCircle = plt.Circle((centroid[0],centroid[1]),radius,fill=False)
      ax.add_patch(centreCircle)
      return ax

def drawNode(nodeCoords, ax, nodeRadius):
      nodeCircle = plt.Circle((nodeCoords[0],nodeCoords[1]), nodeRadius, fc="white", ec="black")
      ax.add_patch(nodeCircle)
      return ax

def drawCentrePolygonPoints(polygonPoints,plt):
      for point in polygonPoints:
            plt.plot(point[0], point[1], "ro")

#Main program
#nodesInfo = setUpNodeDicts([[-260,220],[260,220],[260,-220],[-260,-220]]) #What if they're out of order? Nodes must be sequentially next to their spacewise neighbours.
#nodesInfo = setUpNodeDicts([[-260,220],[260,220],[0,-220]])
#nodesInfo = setUpNodesDicts([])
#nodesInfo = setUpNodeDicts([[-260,-220],[260,-220],[0,220]])
#print(nodesInfo)

def distRatioBeziersAndCentroid(nodeFrom, nodeTo, centroid, radius): #should give t value?
      circleConnectionPoint = findCircleConnectionPointProportion(nodeFrom,nodeTo, radius, centroid)
      nodeFromDist = mag([nodeFrom["coords"][i] - circleConnectionPoint[i] for i in range(2)])
      nodeToDist = mag([nodeTo["coords"][i] - circleConnectionPoint[i] for i in range(2)])
      ratioDivisor = nodeToDist + nodeFromDist
      if nodeFromDist > nodeToDist:
            t = nodeToDist/ ratioDivisor
            print(t)
      else:
            t = nodeFromDist/ ratioDivisor
            print(t)
      #t = 0.7
      print("ccp", circleConnectionPoint)
      controlPoint = [(circleConnectionPoint[i] - nodeFrom["associatedPolygonPoint"][i] + (2*nodeFrom["associatedPolygonPoint"][i]*t) - (nodeFrom["associatedPolygonPoint"][i]*(t**2)) - (nodeTo["associatedPolygonPoint"][i]*(t**2))) / (2*(t - (t**2))) for i in range(2)]
      plt.plot(controlPoint[0], controlPoint[1], "bo")
      print("cp", controlPoint)
      drawCentreCurve([nodeFrom["associatedPolygonPoint"], controlPoint, nodeTo["associatedPolygonPoint"]])
      print("should give control point", [nodeFrom["associatedPolygonPoint"][i] - (2*t*nodeFrom["associatedPolygonPoint"][i]) + (t**2) + (((2-2*t)*t)*controlPoint[i]) + (nodeTo["associatedPolygonPoint"][i]*(t**2)) for i in range(2)])

def tempBeziers(plt):
      points = [[-260,220], [260,220], [0,20]]
      drawCentrePolygonPoints(points, plt)


fig, ax = setUpMatplotCanvas()

def drawEdge(radius, nodeRadius, nodesInfo, ax): #2 edgecase and 1
      print("setup")
      middleConnectorSize = pi #why is this pi???
      for node in nodesInfo:
            ax = drawNode(node["coords"], ax, nodeRadius)
      centroid = findCentroid(nodesInfo)
      #centroid = [-50,-50]
      drawCentreCircle(centroid,radius, ax)
      polygonPoints = calculateCentrePolygonPointsProportion(nodesInfo, centroid, radius, middleConnectorSize)
      drawCentrePolygonPoints(polygonPoints,plt)
      for nodeFromIndex in range(len(nodesInfo)):
            nodeFrom = nodesInfo[nodeFromIndex]
            nodeTo = nodesInfo[(nodeFromIndex + 1) % len(nodesInfo)]
            distRatioBeziersAndCentroid(nodeFrom, nodeTo, centroid, radius)
            #curveVerts = calcCentreCurveVertsProportion(nodeFrom, nodeTo, radius, centroid)
            #drawCentreCurve(curveVerts)
            drawNodeToPolygonLine(nodeFrom["coords"],nodeFrom["associatedPolygonPoint"])

nodesInfo = setUpNodeDicts([[-260,220],[280,90],[260,-220],[-260,-150]])
#nodesInfo = setUpNodeDicts([[-260,220],[260,220],[0,-220]])
radius = 20
nodeRadius = 30 #2 min (tiedto the *2 bit in ) Why is this pi??????
#drawEdge(radius, nodeRadius, nodesInfo, ax)
##nodesInfo = setUpNodeDicts([[260,220],[500,0],[260,-220]])
##drawEdge(radius,nodeRadius,nodesInfo, ax)
###nodesInfo = setUpNodeDicts([[260,220],[260,-220]])
###drawEdge(radius, nodeRadius, nodesInfo, ax)
##nodesInfo = setUpNodeDicts([[-260,220],[130,500], [260,220]])
drawEdge(radius,nodeRadius,nodesInfo, ax)

#tempBeziers(plt)
plt.show()


