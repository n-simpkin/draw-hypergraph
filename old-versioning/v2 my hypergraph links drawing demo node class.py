import turtle as t
from random import randint
from math import sqrt

#Create randomly spaced nodes
##node1Centre = [randint(-270,270) for i in range(2)]
##print(node1Centre)

def drawNode(nodeCoords):
      nodeRadius = 20
      t.teleport(nodeCoords[0],nodeCoords[1])
      t.circle(nodeRadius)
      
def setUpNodeDicts(nodeCoords):
      nodesInfo = []
      properties = ["coords","unitVector","adjustedUnitVector","associatedPolygonPoint"]
      keys = [None for i in range(5)]
      for i,node in enumerate(nodeCoords):
            nodesInfo.append(dict(zip(properties,keys)))
            nodesInfo[i]["coords"] = node
      return nodesInfo

#Find centroid
def findCentroid(nodesInfo):
      centroid = [0,0]
      for node in nodesInfo:
            centroid[0] += node["coords"][0]
            centroid[1] += node["coords"][1]
      centroid[0] /= 3
      centroid[1] /= 3
      return centroid

def drawCentreCircle(centroid,radius):
      #It looks a bit off but that's just how the circumcentre is.
      centroid[1] -= radius #correction for circle drawing from bottom not middle
      t.teleport(centroid[0],centroid[1])
      t.circle(radius)
      centroid[1] += radius
      t.teleport(centroid[0],centroid[1])

#find points of centre polygon, aligned correctly.
def constructCentrePolygon(nodesInfo,centroid,radius):
      polygonPoints = []
      for node in nodesInfo:
            vector = [node["coords"][i] - centroid[i] for i in range(2)] #minus node position from centroid for x and y. There are a few generators that look like this that will come up, they just do it for both the x and y.
            #print("vector of",node["coords"],"is",vector)
            vectorMag = sqrt((vector[0]**2) + (vector[1])**2) # Look out for negatives here
            node["unitVector"] = [vector[i]/vectorMag for i in range(2)]
            #print("unit vector of",node["coords"],"is", node["unitVector"])
            polygonPoint = [(5*radius)*coord for coord in node["unitVector"]] #multiply unit vector by a multiple of the radius so the shape will sit outside the circle
            polygonPoint = [polygonPoint[i] + centroid[i] for i in range(2)] #Adjust the point so it's in relation to the centroid
            node["associatedPolygonPoint"] = polygonPoint
            #print("closestPolygonPoint of",node["coords"],"is", node["associatedPolygonPoint"])
            polygonPoints.append(polygonPoint)
      return polygonPoints
            
def drawCentrePolygon(polygonPoints):
      t.teleport(polygonPoints[0][0],polygonPoints[0][1]) 
      for point in polygonPoints:
            t.goto(point)
      t.goto(polygonPoints[0]) #to complete the shape

def findCircleConnectionPoint(node1, node2, radius, centroid):
##      direction = [None,None]
##      print("unit vector for node",nodes[currentNode],"is",unitVector)
##      print("unit vector for node",nodes[currentNode+1],"is",unitVector)
##      direction[0] = unitVectors[currentNode][0] + unitVectors[currentNode + 1][0]
##      direction[1] = unitVectors[currentNode][1] + unitVectors[currentNode + 1][1]
##      print("direction vector for closest circle point of",currentNode,"is",unitVector)
##      circlePoint = [radius*coord for coord in direction]
##      #Adjust the point so it's in relation to the circumcentre
##      circlePoint[0] = circlePoint[0] + centroid[0]
##      circlePoint[1] = circlePoint[1] + centroid[1]
##      return circlePoint
      direction = [node1["unitVector"][i] + node2["unitVector"][i] for i in range(2)]
      circlePoint = [radius*coord for coord in direction]
      circlePoint = [circlePoint[i] + centroid[i] for i in range(2)]
      print("closest circle point of",node1["coords"], "and", node2["coords"],"is",circlePoint)
      return circlePoint

#def calcLineEquation(nodesInfo):
      

      
#draw the lines
def drawLines(nodesInfo, radius, centroid):
      for nodeIndex in range(len(nodesInfo)-1):
            circleConnectionPoint = findCircleConnectionPoint(nodesInfo[nodeIndex],nodesInfo[nodeIndex+1], radius, centroid)
            t.teleport(circleConnectionPoint[0], circleConnectionPoint[1])


nodesInfo = setUpNodeDicts([[-260,220],[260,220],[0,-220]])
#nodesInfo = setUpNodeDicts([[-260,180],[230,200],[-100,-230]])
radius = 10

print(nodesInfo)
for node in nodesInfo:
      print(node["coords"])
      drawNode(node["coords"])
centroid = findCentroid(nodesInfo)
drawCentreCircle(centroid,radius)
polygonPoints = constructCentrePolygon(nodesInfo, centroid, radius)
print(polygonPoints)
drawCentrePolygon(polygonPoints)
drawLines(nodesInfo, radius, centroid)



