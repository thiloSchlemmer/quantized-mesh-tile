# -*- coding: utf-8 -*-

from __future__ import print_function

from quantized_mesh_tile.editable_terrain import EditableTerrainTile
from quantized_mesh_tile.global_geodetic import GlobalGeodetic

from . import cartesian3d as c3d


def getNextByKeyAndValue(edgeConnections, index, edgeSide):
    """

    :param edgeConnections:
    :param index:
    :param edgeSide:
    :return:
    """
    if index == len(edgeConnections):
        return edgeConnections[index]

    edgeConnectionNext = edgeConnections[index]
    while not edgeConnectionNext.isSide(edgeSide):
        index += 1
        edgeConnectionNext = edgeConnections[index]

    return edgeConnectionNext


def getPreviousByKeyAndValue(edgeConnections, index, edgeSide):
    """

    :param edgeConnections:
    :param index:
    :param edgeSide:
    :return:
    """
    if index == 0:
        return edgeConnections[index]

    edgeConnectionPrevious = edgeConnections[index]
    while not edgeConnectionPrevious.isSide(edgeSide):
        index -= 1
        edgeConnectionPrevious = edgeConnections[index]

    return edgeConnectionPrevious


def getNeighbours(z, x, y):
    return {'west': (z, x - 1, y),
            'north': (z, x, y + 1),
            'south': (z, x, y - 1),
            'east': (z, x + 1, y)}


def getNeighboursSouthEast(z, x, y):
    return {'south': (z, x, y - 1),
            'east': (z, x + 1, y)}


def loadTile(terrainPath, x, y, z):
    """

    :rtype: EditableTerrainTile
    """
    geodetic = GlobalGeodetic(True)
    [minx, miny, maxx, maxy] = geodetic.TileBounds(x, y, z)
    tile = EditableTerrainTile(west=minx, south=miny, east=maxx, north=maxy)
    tile.fromFile(terrainPath, hasLighting=True)
    return tile


class EdgeConnection(object):
    """
    Property class, to store information of points/nodes which
    participating on a tile edge
    """
    BOTH_SIDES = 2
    ONE_SIDE = 1

    def __init__(self, edgeInfo, edgeIndex):
        self.edgeIndex = edgeIndex
        self.edgeInfo = edgeInfo
        self._sideVertices = {}

    def __repr__(self):
        msg = 'E:{0} [{1}] -> ({2})'.format(self.edgeInfo, self.edgeIndex,
                                            self._sideVertices)
        return msg

    def addSide(self, edgeSide, sideVertex):
        """
        Adds the side information to the edge-connection
        :param edgeSide: the participating side of the edge ('w','n','e','s')
        :param sideVertex: the index of the vertex in the list of vertices of the
        given side (tile)
        """
        self._sideVertices[edgeSide] = sideVertex

    def getSideVertex(self, edgeSide):
        """
        Gets the index of the vertex in the list of vertices of the given side (tile)
        :rtype: integer
        :param edgeSide: the participating side of the edge ('w','n','e','s')
        :return: the index of the vertex, based on the given side
        """
        return self._sideVertices[edgeSide]

    @property
    def sideVertices(self):
        return dict(self._sideVertices)

    @property
    def isComplete(self):
        size = len(self._sideVertices.values())
        return EdgeConnection.BOTH_SIDES == size

    @property
    def isBrokenOnCenter(self):
        return self.isBrokenOn(self.edgeInfo)

    @property
    def isBrokenOnNeighbour(self):
        return self.isBrokenOn('c')

    def isBrokenOn(self, edgeSide):
        # type: (str) -> bool
        size = len(self._sideVertices.values())
        return EdgeConnection.ONE_SIDE == size and self.isSide(edgeSide)

    def isSide(self, edgeSide):
        # type: (str) -> bool
        """

        :param edgeSide: the participating side of the edge ('w','n','e','s')
        :return: Returns True, if a vertex of the given side is registered in this
        connection
        """
        return edgeSide in self._sideVertices.keys()


class TileStitcher(object):
    """
        The worker class to stitch terrain files together

        Constructor arguments:

        ''center_tile''
            the the center tile, from which the neighbouring edges are stitched,
            if neighbour tiles are added


        Usage example::
        import os
        from quantized_mesh_tile.tile_stitcher import TileStitcher, loadTile,
                                                        getNeighboursSouthEast

        directory_base_path = '/data/terrain/'
        levels = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

        #walk through level and file hierarchy
        for level in levels:
            directory_path = os.path.join(directory_base_path, str(level))
            terrain_files = []
            for root, dirs, files in os.walk(directory_path, topdown=True):
                for name in files:
                    candidate_path = os.path.join(root, name)
                    if candidate_path.endswith('.terrain'):
                        terrain_files.append(candidate_path)

            for tile_path in terrain_files:

                y = int(os.path.basename(tile_path).split('.')[0])
                x = int(os.path.basename(os.path.dirname(tile_path)))
                print('processing {0} ...'.format(tile_path))
                neighbours = getNeighboursSouthEast(level, x, y)
                center_tile = tile_stitcher.loadTile(tile_path, x, y, level)

                stitcher = TileStitcher(center_tile)
                for n, tile_info in neighbours.items():
                    n_z, n_x, n_y = tile_info

                    neighbour_path = os.path.join(directory_path,
                                    '%s/%s.terrain' % (n_x, n_y))
                    if os.path.exists(neighbour_path):
                        print("\tadding Neighbour {0}...".format(neighbour_path))
                        tile = tile_stitcher.loadTile(neighbour_path, n_x, n_y, level)
                        stitcher.addNeighbour(tile)
                stitcher.stitchTogether()
                stitcher.save()
    """

    def __init__(self, centerTile):
        self._center = centerTile
        self._neighbours = {}

    def _getEdgeConnection(self, neighbourTile):
        centerBbox = self._center.boundingBox
        neighbourBbox = neighbourTile.boundingBox

        if centerBbox['west'] == neighbourBbox['east']:
            return 'w'
        if centerBbox['east'] == neighbourBbox['west']:
            return 'e'
        if centerBbox['north'] == neighbourBbox['south']:
            return 'n'
        if centerBbox['south'] == neighbourBbox['north']:
            return 's'
        return None

    def _getEdgeIndices(self, neighbourTile):
        centerBbox = self._center.boundingBox
        neighbourBbox = neighbourTile.boundingBox
        centerVertices = []
        neighbourVertices = []
        if centerBbox['west'] == neighbourBbox['east']:
            centerVertices = self._center.getEdgeIndices('w')
            neighbourVertices = neighbourTile.getEdgeIndices('e')
        if centerBbox['east'] == neighbourBbox['west']:
            centerVertices = self._center.getEdgeIndices('e')
            neighbourVertices = neighbourTile.getEdgeIndices('w')
        if centerBbox['north'] == neighbourBbox['south']:
            centerVertices = self._center.getEdgeIndices('n')
            neighbourVertices = neighbourTile.getEdgeIndices('s')
        if centerBbox['south'] == neighbourBbox['north']:
            centerVertices = self._center.getEdgeIndices('s')
            neighbourVertices = neighbourTile.getEdgeIndices('n')

        return centerVertices, neighbourVertices

    def _findEdgeConnections(self):

        edgeConnections = []
        for edgeInfo, neighbourTile in self._neighbours.items():
            singleEdgeVertices = {}
            edgeIndex = 0  # assume south>|<north
            if edgeInfo in ['w', 'e']:  # assume west>|<east
                edgeIndex = 1

            centerIndices, neighbourIndices = self._getEdgeIndices(neighbourTile)
            for centerIndex in centerIndices:
                cUV = (self._center.u[centerIndex], self._center.v[centerIndex])
                cKey = '{}_{:05}'.format(edgeInfo, cUV[edgeIndex])
                edgeConnection = EdgeConnection(edgeInfo, cUV[edgeIndex])
                edgeConnection.addSide('c', centerIndex)
                singleEdgeVertices[cKey] = edgeConnection
            for neighbourIndex in neighbourIndices:
                nUV = (neighbourTile.u[neighbourIndex],
                       neighbourTile.v[neighbourIndex])
                nKey = '{}_{:05}'.format(edgeInfo, nUV[edgeIndex])

                if nKey in singleEdgeVertices.keys():
                    singleEdgeVertices[nKey].addSide(edgeInfo, neighbourIndex)
                else:
                    edgeConnection = EdgeConnection(edgeInfo, nUV[edgeIndex])
                    edgeConnection.addSide(edgeInfo, neighbourIndex)
                    singleEdgeVertices[nKey] = edgeConnection

            singleEdgeConnections = sorted(singleEdgeVertices.values(),
                                           key=lambda x: x.edgeIndex,
                                           reverse=False)

            edgeConnections.append(singleEdgeConnections)

        return edgeConnections

    def _stitchEdges(self, edgeConnections):

        for edge in edgeConnections:
            for index in range(len(edge)):
                edgeConnection = edge[index]

                edgeInfo = edgeConnection.edgeInfo
                neighbour = self._neighbours[edgeInfo]
                if edgeConnection.isComplete:
                    self._updateHeightToEven(edgeConnection)
                elif edgeConnection.isBrokenOnNeighbour:
                    previousVertex = self._getPreviousVertex(index, edge, edgeInfo)
                    nextVertex = self._getNextVertex(index, edge, edgeInfo)

                    splittingCoordinate = self._center.getCoordinateAt(
                        edgeConnection.getSideVertex('c'))
                    newVertex = neighbour.findAndSplitTriangle(previousVertex,
                                                               nextVertex,
                                                               splittingCoordinate)
                    edgeConnection.addSide(edgeInfo, newVertex)
                else:
                    # wenn vertex nur in n, dann triangle in c
                    # von c-vertex-1 und c-vertex+1 splitten
                    previousVertex = self._getPreviousVertex(index, edge, 'c')
                    nextVertex = self._getNextVertex(index, edge, 'c')

                    splittingCoordinate = neighbour.getCoordinateAt(
                        edgeConnection.getSideVertex(edgeInfo))
                    newVertex = self._center.findAndSplitTriangle(previousVertex,
                                                                  nextVertex,
                                                                  splittingCoordinate)

                    edgeConnection.addSide('c', newVertex)

    def _harmonizeNormals(self, edgeConnections):
        center = self._center
        for edge in edgeConnections:
            for edgeConnection in edge:
                centerVertexIndex = edgeConnection.getSideVertex('c')

                sideVertex = edgeConnection.getSideVertex(edgeConnection.edgeInfo)
                neighbourVertexIndices = {edgeConnection.edgeInfo: sideVertex}

                centerTriangles = center.findAllTrianglesOf(centerVertexIndex)
                normals = center.calculateWeightedNormalsFor(centerTriangles)

                for neighbourInfo, vertexIndex in neighbourVertexIndices.items():
                    neighbourTile = self._neighbours[neighbourInfo]
                    triangles = neighbourTile.findAllTrianglesOf(vertexIndex)
                    normals.extend(neighbourTile.calculateWeightedNormalsFor(
                        triangles))

                normalVertex = [0, 0, 0]
                for weightedNormal in normals:
                    normalVertex = c3d.add(normalVertex, weightedNormal)

                normalVertex = c3d.normalize(normalVertex)
                center.setNormal(centerVertexIndex, normalVertex)
                for neighbourInfo, vertexIndex in neighbourVertexIndices.items():
                    neighbourTile = self._neighbours[neighbourInfo]
                    neighbourTile.setNormal(vertexIndex, normalVertex)

    def _buildNormals(self, edgeConnections):
        center = self._center
        center.rebuildH()
        for n in self._neighbours.values():
            n.rebuildH()

        for edge in edgeConnections:
            for edgeConnection in edge:
                centerVertexIndex = edgeConnection.getSideVertex('c')

                sideVertex = edgeConnection.getSideVertex(edgeConnection.edgeInfo)
                neighbourVertexIndices = {edgeConnection.edgeInfo: sideVertex}

                centerTriangles = center.findAllTrianglesOf(centerVertexIndex)
                normals = center.calculateWeightedNormalsFor(centerTriangles)

                for neighbourInfo, vertexIndex in neighbourVertexIndices.items():
                    neighbourTile = self._neighbours[neighbourInfo]
                    triangles = neighbourTile.findAllTrianglesOf(vertexIndex)
                    normals.extend(neighbourTile.calculateWeightedNormalsFor(
                        triangles))

                normalVertex = [0, 0, 0]
                for weightedNormal in normals:
                    normalVertex = c3d.add(normalVertex, weightedNormal)

                normalVertex = c3d.normalize(normalVertex)
                center.setNormal(centerVertexIndex, normalVertex)
                for neighbourInfo, vertexIndex in neighbourVertexIndices.items():
                    neighbourTile = self._neighbours[neighbourInfo]
                    neighbourTile.setNormal(vertexIndex, normalVertex)

    @staticmethod
    def _getNextVertex(index, edgeConnections, edgeSide):
        edgeConnectionNext = getNextByKeyAndValue(edgeConnections,
                                                  index,
                                                  edgeSide)
        vertexNext = edgeConnectionNext.getSideVertex(edgeSide)
        return vertexNext

    @staticmethod
    def _getPreviousVertex(index, edgeConnections, edgeSide):
        edgeConnectionPrevious = getPreviousByKeyAndValue(edgeConnections,
                                                          index,
                                                          edgeSide)
        vertexPrevious = edgeConnectionPrevious.getSideVertex(edgeSide)
        return vertexPrevious

    def _updateHeightToEven(self, edgeConnection):
        centerVertexIndex = edgeConnection.getSideVertex('c')

        vertexIndices = edgeConnection.sideVertices

        vertexHeights = []
        for edgeInfo, vertexIndex in vertexIndices.items():
            if edgeInfo is 'c':
                vertexHeights.append(self._center.getHeightAt(vertexIndex))
            else:
                neighbour = self._neighbours[edgeInfo]
                vertexHeights.append(neighbour.getHeightAt(vertexIndex))

        height = sum(vertexHeights) / len(vertexHeights)
        for edgeInfo in vertexIndices:
            if edgeInfo is 'c':
                self._center.setHeight(centerVertexIndex, height)
            else:
                neighbour = self._neighbours[edgeInfo]
                neighbour.setHeight(vertexIndices[edgeInfo], height)

    def addNeighbour(self, neighbourTile):
        edgeConnection = self._getEdgeConnection(neighbourTile)
        self._neighbours[edgeConnection] = neighbourTile

    def harmonizeNormals(self):
        edgeConnections = self._findEdgeConnections()
        self._harmonizeNormals(edgeConnections)

    def stitchTogether(self):
        edgeConnections = self._findEdgeConnections()
        self._stitchEdges(edgeConnections)
        self._buildNormals(edgeConnections)

    def save(self):
        self._center.save()
        for edgeInfo, neighbourTile in self._neighbours.items():
            neighbourTile.save()

    def saveTo(self, dir_path):
        self._center.saveTo(dir_path)
        for edgeInfo, neighbourTile in self._neighbours.items():
            neighbourTile.saveTo(dir_path)
