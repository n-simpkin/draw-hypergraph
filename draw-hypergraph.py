from math import atan, cos, degrees, pi, radians, sin, sqrt

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np  # rewrite all wiht numpy arrays, eliminates all my clunky two generators, much more semantic.
from matplotlib.path import Path

"""
TODOS
Use numpy for higher readability from clunky generators
Factor out magnitude and unit vector
Consider structure
Hypergraph class for holding the information?

"""


# SETUP


def setUpMatplotCanvas():
    fig, ax = plt.subplots()
    plt.xlim(-300, 300)
    plt.ylim(-300, 300)
    ax.set_aspect(1)
    return fig, ax


def setUpNodeDicts(nodeCoords):
    nodesInfo = []
    keys = ["coords", "unitVector", "associatedPolygonPoint", "t"]
    vals = [None for i in range(4)]

    for i, node in enumerate(np.array(nodeCoords)):
        nodesInfo.append(dict(zip(keys, vals)))
        nodesInfo[i]["coords"] = node
    return nodesInfo


# CALCULATION TOOLS


def magnitude(coords):
    return sqrt((coords[0] ** 2) + (coords[1] ** 2))


def makeUnitVector(vector):
    return np.divide(vector, magnitude(vector))


# CALCULATE


def findCentroid(nodesInfo):
    # Take mean average of all coordinates to find the centroid
    centroid = np.array([0, 0])

    for node in nodesInfo:
        centroid = np.add(node["coords"], centroid)
    centroid = centroid / 3

    return centroid


def calculatePolygonPoints(
    nodesInfo, centroid, radius, polygonPointDistance
):  # find points of centre polygon, aligned correctly.
    polygonPoints = []

    for node in nodesInfo:
        vector = np.subtract(node["coords"], centroid)
        node["unitVector"] = makeUnitVector(
            vector
        )  # Get direction from the centroid to the node
        polygonPoint = np.multiply(
            vector, polygonPointDistance
        )  # multiply unit vector by a multiple of the radius so the shape will sit outside the circle
        polygonPoint = np.add(
            polygonPoint, centroid
        )  # translate the point so its centred on the centroid rather than the origin
        node["associatedPolygonPoint"] = polygonPoint
        polygonPoints.append(polygonPoint)

    return polygonPoints


def calculateCircleConnectionPoint(node1, node2, radius, centroid):  # Should be fixed
    resultantVector = np.add(
        node1["unitVector"], node2["unitVector"]
    )  # Doesn't produce a unit vector so the point isn't actually radius distance away.
    direction = makeUnitVector(resultantVector)
    circlePoint = np.multiply(
        direction, radius
    )  # The multiply the direction (of length one) by the circle radius to find the point on the circle needed.
    circlePoint = np.add(
        circlePoint, centroid
    )  # All this will have been performed with the centre as 0,0, so translate the point so its in relation to the centorid
    plt.plot(circlePoint[0], circlePoint[1], "yo")
    return circlePoint


def calculateBezierControlPoint(
    nodeFrom, nodeTo, centroid, radius
):  # should give t value?
    circleConnectionPoint = calculateCircleConnectionPoint(
        nodeFrom, nodeTo, radius, centroid
    )
    nodeFromDist = magnitude(np.subtract(nodeFrom["coords"], circleConnectionPoint))
    nodeToDist = magnitude(np.subtract(nodeTo["coords"], circleConnectionPoint))
    ratioDivisor = nodeToDist + nodeFromDist

    # To get the right t value I need the smaller value divided by the sum, so just chekcing for that.
    if nodeFromDist > nodeToDist:
        t = nodeToDist / ratioDivisor
    else:
        t = nodeFromDist / ratioDivisor
    nodeFrom["t"] = t

    # Calculates where to put the control point if I want the curve to start at nodeFrom polygon point, end at nodeTo polygon point, and have the turning point at the circle connection point.
    # Theory for this is written in full in my first notebook and mainly sourced from https://pomax.github.io/bezierinfo/
    # Bezier equation for a quadratic is (where x(t) is a function) x(t) = X1(1-t)^2 + (X2)(2)(1-t)(t) + (X3)(t^2).
    # Rearranged for X3 - The control point - this is X2 = (x(t) - X1(1-t)^2 - (X3)(t^2))/ (2)(1-t)(t)
    return [
        (
            circleConnectionPoint[i]
            - nodeFrom["associatedPolygonPoint"][i]
            + (2 * nodeFrom["associatedPolygonPoint"][i] * t)
            - (nodeFrom["associatedPolygonPoint"][i] * (t**2))
            - (nodeTo["associatedPolygonPoint"][i] * (t**2))
        )
        / (2 * (t - (t**2)))
        for i in range(2)
    ]


# def calculateBezierRatios(nodeFrom, nodeTo, centroid):

#     return [nodeFrom["associatedPolygonPoints"], centroid, nodeTo["associatedPolygonPoints"]]


def calculateRatioBetweenNodesAndCentroid(nodeFromCoords, nodeToCoords, centroid):
    nodeFromDist = magnitude(np.subtract(nodeFromCoords, centroid))
    nodeToDist = magnitude(np.subtract(nodeToCoords, centroid))
    ratioDivisor = nodeToDist + nodeFromDist

    # To get the right t value I need the smaller value divided by the sum, so just chekcing for that.
    if nodeFromDist > nodeToDist:
        nodesToCentroidDistanceRatio = nodeToDist / ratioDivisor
    else:
        nodesToCentroidDistanceRatio = nodeFromDist / ratioDivisor

    return nodesToCentroidDistanceRatio


def calculateRatioPointBetweenNodes(
    nodeFromCoords, nodeToCoords, nodesToCentroidDistanceRatio
):
    nodeToNodeVector = np.subtract(nodeToCoords, nodeFromCoords)  # This is right
    nodeToNodeDist = magnitude(nodeToNodeVector)
    nodeToNodeDirection = makeUnitVector(nodeToNodeVector)

    ratioPoint = nodeToNodeDirection * (nodeToNodeDist * nodesToCentroidDistanceRatio)
    ratioPoint = nodeFromCoords + ratioPoint
    return ratioPoint


def RationalBezier(t, weights, ratios):
    # This is basically word for word from https://pomax.github.io/bezierinfo/#weightcontrol
    # w = [-260, 220], [260, 220], [0, -220]
    # r = [1, 2, 0.8]
    t2 = t**2
    mt = 1 - t
    mt2 = mt**2
    f = [ratios[0] * mt2, 2 * ratios[1] * mt * t, ratios[2] * t2]
    basis = f[0] + f[1] + f[2]
    return (
        f[0] * weights[0] + f[1] * weights[1] + f[2] * weights[2]
    ) / basis  # Gives the coordinates of the given curve at the given t value.


# DRAW


def drawCentreCurve(curveVerts):
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    bezier1 = patches.PathPatch(Path(curveVerts, codes), fc="none")
    ax.add_patch(bezier1)


def drawNodeToPolygonLine(nodeCoords, polygonPointCoords):
    plt.plot(
        (nodeCoords[0], polygonPointCoords[0]),
        (nodeCoords[1], polygonPointCoords[1]),
        color="black",
    )


def drawDashedLine(pointFrom, pointTo):
    plt.plot((pointFrom[0], pointTo[0]), (pointFrom[1], pointTo[1]), "--")


def drawCircle(coords, radius, ax):
    circle = plt.Circle((coords[0], coords[1]), radius, fc="white", ec="black")
    ax.add_patch(circle)
    return ax


def drawPoints(points, plt):
    for point in points:
        plt.plot(point[0], point[1], "ro")


fig, ax = setUpMatplotCanvas()


def drawEdge(radius, nodeRadius, nodesInfo, ax):  # 2 edgecase and 1
    polygonPointDistance = 1

    for node in nodesInfo:
        ax = drawCircle(node["coords"], nodeRadius, ax)
    centroid = findCentroid(nodesInfo)
    drawCircle(centroid, radius, ax)
    polygonPoints = calculatePolygonPoints(
        nodesInfo, centroid, radius, polygonPointDistance
    )
    drawPoints(polygonPoints, plt)

    for nodeFromIndex in range(len(nodesInfo)):
        nodeFrom = nodesInfo[nodeFromIndex]
        nodeTo = nodesInfo[(nodeFromIndex + 1) % len(nodesInfo)]
        # controlPoint = calculateBezierControlPoint(nodeFrom, nodeTo, centroid, radius)
        controlPoint = centroid
        # plt.plot(controlPoint[0], controlPoint[1], "bo")
        nodesToCentroidDistanceRatio = calculateRatioBetweenNodesAndCentroid(
            nodeFrom["coords"], nodeTo["coords"], centroid
        )
        # calculateRatioPointBetweenNodes(nodeFrom["coords"], nodeTo["coords"], nodesToCentroidDistanceRatio)
        ratioPoint = calculateRatioPointBetweenNodes(
            nodeFrom["coords"], nodeTo["coords"], nodesToCentroidDistanceRatio
        )

        drawCentreCurve(
            [
                nodeFrom["associatedPolygonPoint"],
                controlPoint,
                nodeTo["associatedPolygonPoint"],
            ]
        )
        drawNodeToPolygonLine(nodeFrom["coords"], nodeFrom["associatedPolygonPoint"])
        drawDashedLine(nodeFrom["coords"], nodeTo["coords"])
        drawDashedLine(ratioPoint, centroid)
        drawDashedLine(nodeFrom["coords"], centroid)
        drawPoints([ratioPoint], plt)


nodesInfo = setUpNodeDicts([[-260, 220], [90, 90], [260, -220], [-260, -150]])
# nodesInfo = setUpNodeDicts([[-260, 220], [260, 220], [0, -220]])
# nodesInfo = setUpNodeDicts([[260, 220], [500, 0], [260, -220]])
###nodesInfo = setUpNodeDicts([[260,220],[260,-220]])
# nodesInfo = setUpNodeDicts([[-260, 220], [130, 500], [260, 220]])

radius = 25
nodeRadius = 30

drawEdge(radius, nodeRadius, nodesInfo, ax)
plt.show()
