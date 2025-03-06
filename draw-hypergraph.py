from math import atan, cos, degrees, pi, radians, sin, sqrt

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np  # rewrite all wiht numpy arrays, eliminates all my clunky two generators, much more semantic.
from matplotlib.path import Path
from matplotlib.widgets import Slider

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
    plt.subplots_adjust(bottom=0.2)
    plt.xlim(-300, 300)
    plt.ylim(-300, 300)
    ax.set_aspect(1)
    return fig, ax


def setUpNodeDicts(nodeCoords):
    nodesInfo = []
    keys = ["coords", "unitVector", "associatedPolygonPoint"]
    vals = [None for i in range(3)]

    for i, node in enumerate(np.array(nodeCoords)):
        nodesInfo.append(dict(zip(keys, vals)))
        nodesInfo[i]["coords"] = node
    return nodesInfo


def setUpBeziersDicts(numOfCurves):
    beziersInfo = []
    keys = ["weights", "ratios"]
    vals = [None for i in range(2)]

    for i in range(numOfCurves):
        beziersInfo.append(dict(zip(keys, vals)))
    return beziersInfo


# CALCULATION TOOLS


def magnitude(coords):
    return sqrt((coords[0] ** 2) + (coords[1] ** 2))


def makeUnitVector(vector):
    return np.divide(vector, magnitude(vector))


# CALCULATE


def findCentroid(points):
    # Take mean average of all coordinates to find the centroid
    centroid = np.array([0, 0])

    for point in points:
        centroid = np.add(point, centroid)
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


def calcRationalBezierPoint(t, weights, ratios):
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


def calculateBezierPlotPointsBySegments(
    bezierInfo, segmentCount
):  # Simple approximate curve drawing. Tempoary unitl I can be bothered to make a better one.
    step = 1 / segmentCount
    # coords = [bezierInfo["weights"][0]]  # Starts coords at start point
    coords = []
    for i in range(segmentCount + 1):
        t = step * i
        coords.append(
            calcRationalBezierPoint(t, bezierInfo["weights"], bezierInfo["ratios"])
        )
    return coords


def calculateTP(
    nodeFromCoords, nodeToCoords, centroid, bezierInfo
):  # Using centroid is a bit wrong as should be circle connecion point, but I don't have that anymore and I'm a bit lazy.
    nodeFromDist = magnitude(np.subtract(nodeFromCoords, centroid))
    nodeToDist = magnitude(np.subtract(nodeToCoords, centroid))
    ratioDivisor = nodeToDist + nodeFromDist

    if nodeFromDist > nodeToDist:
        t = nodeToDist / ratioDivisor
    else:
        t = nodeFromDist / ratioDivisor

    coords = calcRationalBezierPoint(t, bezierInfo["weights"], bezierInfo["ratios"])

    return coords


# DRAW


# def addSlider():
#     axRatio = fig.add_axes([0.25, 0.1, 0.65, 0.03])
#     ratioSlider = Slider(
#         ax=axRatio,
#         label="Control point ratio",
#         valmin=0,
#         valmax=2,
#         valinit=1,
#     )

# def update(val):
#     ratio = ratioSlider.val


def drawNodeToPolygonLine(nodeCoords, polygonPointCoords):
    plt.plot(
        (nodeCoords[0], polygonPointCoords[0]),
        (nodeCoords[1], polygonPointCoords[1]),
        color="black",
    )


def drawDashedLine(pointFrom, pointTo):
    plt.plot((pointFrom[0], pointTo[0]), (pointFrom[1], pointTo[1]), "--")


def drawCircle(coords, radius):
    circle = plt.Circle((coords[0], coords[1]), radius, fc="white", ec="black")
    ax.add_patch(circle)


def drawPoints(points, colour="ro"):
    for point in points:
        (pointObj,) = ax.plot(point[0], point[1], colour)
    return pointObj


def drawBezierBySegments(coords):
    Xs = [coordPair[0] for coordPair in coords]
    Ys = [coordPair[1] for coordPair in coords]
    (line,) = ax.plot(Xs, Ys)
    return line


fig, ax = setUpMatplotCanvas()


def drawEdge(
    radius, nodeRadius, nodesList, polygonPointDistance, centre
):  # 2 edgecase and 1
    nodesInfo = setUpNodeDicts(nodesList)
    beziersInfo = setUpBeziersDicts(len(nodesInfo))
    centroid = centre  # Update so everything is centre, cba rn.

    for node in nodesInfo:
        drawCircle(node["coords"], nodeRadius)

    drawCircle(centroid, radius)
    polygonPoints = calculatePolygonPoints(
        nodesInfo, centroid, radius, polygonPointDistance
    )
    drawPoints(polygonPoints)

    linesManipulatable = []
    turningPointsManipulatable = []

    for nodeFromIndex in range(len(nodesInfo)):
        nodeFrom = nodesInfo[nodeFromIndex]
        nodeTo = nodesInfo[(nodeFromIndex + 1) % len(nodesInfo)]
        beziersInfo[nodeFromIndex]["weights"] = (
            nodeFrom["associatedPolygonPoint"],
            centroid,
            nodeTo["associatedPolygonPoint"],
        )
        beziersInfo[nodeFromIndex]["ratios"] = [1, 1, 1]  # Initialise as normal bezier

        nodesToCentroidDistanceRatio = calculateRatioBetweenNodesAndCentroid(
            nodeFrom["coords"], nodeTo["coords"], centroid
        )
        calculateRatioPointBetweenNodes(
            nodeFrom["coords"], nodeTo["coords"], nodesToCentroidDistanceRatio
        )
        ratioPoint = calculateRatioPointBetweenNodes(
            nodeFrom["coords"], nodeTo["coords"], nodesToCentroidDistanceRatio
        )
        turningPoint = calculateTP(
            nodeFrom["coords"], nodeTo["coords"], centroid, beziersInfo[nodeFromIndex]
        )

        bezierCoords = calculateBezierPlotPointsBySegments(
            beziersInfo[nodeFromIndex], 40
        )
        linesManipulatable.append(drawBezierBySegments(bezierCoords))

        drawNodeToPolygonLine(nodeFrom["coords"], nodeFrom["associatedPolygonPoint"])
        drawDashedLine(nodeFrom["coords"], nodeTo["coords"])
        drawDashedLine(ratioPoint, centroid)
        drawDashedLine(nodeFrom["coords"], centroid)

        drawPoints([ratioPoint])
        turningPointsManipulatable.append(drawPoints([turningPoint], "go"))
    return linesManipulatable, beziersInfo, turningPointsManipulatable


# nodesList = [[-260, 220], [90, 90], [260, -220], [-260, -150]]
nodesList = [[-260, 220], [260, 220], [0, -220]]
# nodesList = [[260, 220], [100, 0], [260, -220]]
# nodesList = [[260,220],[260,-220]]
# nodesList = [[-260, 220], [130, 500], [260, 220]]

centroid = findCentroid(nodesList)
centre = centroid

radius = 25
nodeRadius = 30
polygonPointDistance = 1

lines, beziersInfo, turningPoints = drawEdge(
    radius, nodeRadius, nodesList, polygonPointDistance, centre
)
lastCurveTemp = 0


axRatio = fig.add_axes([0.25, 0.1, 0.65, 0.03])
ratioSlider = Slider(
    ax=axRatio,
    label="Control point ratio",
    valmin=0.00001,
    valmax=2,
    valinit=1,
)
axRatio2 = fig.add_axes([0.25, 0.05, 0.65, 0.03])
ratioSlider2 = Slider(
    ax=axRatio2,
    label="Other points ratio",
    valmin=0.00001,
    valmax=2,
    valinit=1,
)


def update(val):
    ratioControlPoint = ratioSlider.val
    ratioOther = ratioSlider2.val
    for i, line in enumerate(lines):
        beziersInfo[i]["ratios"] = [
            ratioOther,
            ratioControlPoint,
            ratioOther,
        ]
        coords = calculateBezierPlotPointsBySegments(beziersInfo[i], 40)
        Xs = [coordPair[0] for coordPair in coords]
        Ys = [coordPair[1] for coordPair in coords]
        line.set_xdata(Xs)
        line.set_ydata(Ys)

    for i, tpObj in enumerate(turningPoints):
        tp = calculateTP(
            nodesList[i], nodesList[(i + 1) % len(nodesList)], centre, beziersInfo[i]
        )
        tpObj.set_xdata([tp[0]])
        tpObj.set_ydata([tp[1]])


ratioSlider.on_changed(update)
ratioSlider2.on_changed(update)

plt.show()
