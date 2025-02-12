import turtle as t
from random import randint
from math import sqrt

#Create randomly spaced nodes
##node1Centre = [randint(-270,270) for i in range(2)]
##print(node1Centre)

#Set up nodes
#canvas -300 to 300 x and  -270 to 270 y
def setUpNodes():
      nodes = [[-260,220],[260,220],[0,-220]]
      t.speed(10)
      nodeRadius = 20
      #Draw nodes
      for node in nodes:
            t.teleport(node[0],node[1])
            t.circle(nodeRadius)
      return nodes

nodes = setUpNodes()
#Find centroid
centroid = [0,0]
for node in nodes:
      centroid[0] += node[0]
      centroid[1] += node[1]
centroid[0] /= 3
centroid[1] /= 3

#draw centre little circle
radius = 10
centroid[1] -= radius #correction for circle drawing from bottom not middle
t.teleport(centroid[0],centroid[1])
t.circle(radius)
centroid[1] += radius
t.teleport(centroid[0],centroid[1])
#It looks a bit off but that's just how the circumcentre is.

#find points of centre polygon, aligned correctly.
vector = [0,0]
unitVector = [0,0]
centreShapeVerticies = []
unitVectors = [0 for node in nodes]
print(unitVectors)
for node in nodes:
      #Calculate unit vector
      vector[0] = node[0] - centroid[0]
      vector[1] = node[1] - centroid[1]
      vectorMag = sqrt((vector[0]**2) + (vector[1])**2) # Look out for negatives here
      unitVector[0] = vector[0] / vectorMag
      unitVector[1] = vector[1] / vectorMag
      print("unit vector for node",node,"is",unitVector)
      unitVectors.append(unitVector)
      print("unit vectors:", unitVectors)
      #multiply unit vector by a multiple of the radius so the shape will sit outside the circle
      point = [(5*radius)*coord for coord in unitVector] #assumes non 0
      #Adjust the point so it's in relation to the circumcentre
      point[0] = point[0] + centroid[0]
      point[1] = point[1] + centroid[1]
      centreShapeVerticies.append(point)
      
print(centroid)
print(centreShapeVerticies)
print(unitVectors)

#draw centre polygon
t.teleport(centreShapeVerticies[0][0],centreShapeVerticies[0][1]) 
for point in centreShapeVerticies:
      t.goto(point)
t.goto(centreShapeVerticies[0]) #to complete the shape

def findCircleConnectingPoint(nodes,currentNode):
      direction = [None,None]
      print("unit vector for node",nodes[currentNode],"is",unitVector)
      print("unit vector for node",nodes[currentNode+1],"is",unitVector)
      direction[0] = unitVectors[currentNode][0] + unitVectors[currentNode + 1][0]
      direction[1] = unitVectors[currentNode][1] + unitVectors[currentNode + 1][1]
      print("direction vector for closest circle point of",currentNode,"is",unitVector)
      circlePoint = [radius*coord for coord in direction]
      #Adjust the point so it's in relation to the circumcentre
      circlePoint[0] = circlePoint[0] + centroid[0]
      circlePoint[1] = circlePoint[1] + centroid[1]
      return circlePoint
      
#draw the lines
for node in range(0,len(nodes)-2):
      pointTo = findCircleConnectingPoint(nodes,node)
      print("closest circle point", pointTo)
      t.teleport(pointTo[0],pointTo[1])


##class nodesAndLinkGraph():
##
##      def ___init___(self, nodes):
##            self.nodes = nodes
##            self.centroid = None





