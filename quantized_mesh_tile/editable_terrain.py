# -*- coding: utf-8 -*-
import numpy
import os

import numpy as np
from future.utils import old_div

from quantized_mesh_tile import TerrainTile
from quantized_mesh_tile.llh_ecef import LLH2ECEF
from quantized_mesh_tile.terrain import lerp
from quantized_mesh_tile.utils import triangleArea
from . import cartesian3d as c3d

null_normal = [0, 0, 0]
np_null_normal = np.array([0, 0, 0])


class EditableTerrainTile(TerrainTile):
    """
    Helper class to edit a given terrain tile on the edges.
    Changes are possible on specified Heights, Triangles and Normals.
    """

    def __init__(self, *args, **kwargs):
        super(EditableTerrainTile, self).__init__(*args, **kwargs)
        self.isIndexDirty = False
        self._changedHeights = []
        self._filePath = None
        self._gzipped = False
        self.name = None

    def setName(self, name):
        self.name = name

    def getEdgeIndices(self, edge):
        """
        Returns the indices of vertices on the defined edge.
        :param edge: the edge of the tile. Specified by the flags 'w','n','e','s'
                    for west-, north-, east- and  south-edge
        :return: Array of integers
        """
        if 'w' == edge:
            edgeValue = TerrainTile.MIN
            searchArray = self.u
        elif 'e' == edge:
            edgeValue = TerrainTile.MAX
            searchArray = self.u
        elif 'n' == edge:
            edgeValue = TerrainTile.MAX
            searchArray = self.v
        elif 's' == edge:
            edgeValue = TerrainTile.MIN
            searchArray = self.v
        indices = [i for i, x in enumerate(searchArray) if x == edgeValue]

        if len(indices) == 0:
            raise Exception("No edge vertices found for edge: {}".format(edge))
        return indices

    def getEdgeCoordinates(self, edge):
        # type: (string) -> Array
        """
        Returns an array of coordinate-tupels for all vertices on the specified edge
        :return: Array of 3d coordinates (3-tupel of float; x,y,z)
        :param edge: the edge of the tile. Specified by the flags 'w','n','e','s'
                    for west, north, east and  south
        """
        coordinates = self.getVerticesCoordinates()
        return [coordinates[v] for v in self.getEdgeIndices(edge)]

    @property
    def boundingBox(self):
        """
        Returns the bounding box of the tile in WGS84 degree
        :return: Dictionary of floats for boundingbox of the tile with
                keys ['west','east','north', 'south']
        """
        return {'west': self._west,
                'east': self._east,
                'north': self._north,
                'south': self._south}

    def setNormal(self, index, normal):
        # type: (int, Array) -> void
        """
        Sets the normal vector for the vertex which is specified with the index,
        changing a normal vector causes a rebuild of all indices
        :param index: index of the vertex
        :param normal: the normal vector in the form [x,y,z]
        """
        self.isIndexDirty = True
        self.vLight[index] = normal

    def setHeight(self, index, height):
        # type: (int, float) -> void
        """
        Sets the height for the vertex  which is specified with the index,
        if the new height is out of the tile height range (min/max height),
        all heights are requantized during save-process
        :param index: index of the vertex
        :param height:the height at the specified vertex index
        """
        isHeightDirty = height < self.header['minimumHeight'] or self.header[
            'maximumHeight'] < height

        if isHeightDirty or self._changedHeights:
            if not self._changedHeights:
                self._changedHeights = [self._dequantizeHeight(x) for x in self.h]
            self._changedHeights[index] = height
        else:
            self.h[index] = self._quantizeHeight(height)

    def getHeightAt(self, index):
        # type: (int) -> float
        """
        Returns the dequantized height at the specified vertex index
        :param index: index of the vertex
        :return: the height at the specified vertex
        """
        height = self._dequantizeHeight(self.h[index])
        return height

    def getCoordinateAt(self, index):
        """
        Returns the dequantized coordinate at the specified vertex index
        :param index: index of the vertex
        :return: the wgs84 coordinate in the form [longitude, latitude, height]
        """
        return self._UVH2LLH(index)

    def findTriangleOf(self, indexOfPreviousVertex, indexOfNextVertex):
        """
        Returns the triangle where both given vertices are member of.
        :param indexOfPreviousVertex: index of the first (previous) vertex
        :param indexOfNextVertex: index of the second (next) vertex
        :return: the index of the queried triangle or None
        """
        trianglesOfPrevious = self.findAllTrianglesOf(indexOfPreviousVertex)
        trianglesOfNext = self.findAllTrianglesOf(indexOfNextVertex)

        indices = list(trianglesOfPrevious & trianglesOfNext)
        if indices:
            return indices[0]

        return None

    def findAllTrianglesOf(self, indexOfVertex):
        """
        Searches for all triangles of the specified vertex index
        and returns the list of triangles
        :param indexOfVertex: the vertex index
        :return: Array of triangle indices
        """
        return set([int(i / 3) for i, v in enumerate(self.indices) if v == indexOfVertex])

    def getTriangleAt(self, index):
        """
        Returns the triangle for the specified index
        :param index: the index of the triangle
        :return: array of vertex indeces
        """
        offset = int(index * 3)

        vi1 = self.indices[offset]
        vi2 = self.indices[offset + 1]
        vi3 = self.indices[offset + 2]
        return vi1, vi2, vi3

    def calculateWeightedNormalsFor(self, triangleIndices):
        """
        Calculates normal vectors for the specified triangles, this
        normal vectors are not normalized and multiplicated with the
        area of participating triangleIndex
        :rtype: Array
        :param triangles:
        :return: Array of not normalized vectors [float,float,float]
        """
        weightedNormals = []
        for triangleIndex in triangleIndices:
            i0, i1, i2 = self.getTriangleAt(triangleIndex)
            llh0 = self._UVH2LLH(i0)
            llh1 = self._UVH2LLH(i1)
            llh2 = self._UVH2LLH(i2)
            v0 = LLH2ECEF(llh0[0], llh0[1], llh0[2])
            v1 = LLH2ECEF(llh1[0], llh1[1], llh1[2])
            v2 = LLH2ECEF(llh2[0], llh2[1], llh2[2])

            normal = np.cross(c3d.subtract(v1, v0), c3d.subtract(v2, v0))
            area = triangleArea(v0, v1)
            weightedNormals.append(normal * area)

        return weightedNormals

    def toFile(self, filePath, gzipped=False):
        if self.isIndexDirty:
            self._rebuildIndices()

        super(EditableTerrainTile, self).toFile(filePath, gzipped)

    def fromFile(self, filePath, hasLighting=False, hasWatermask=False, gzipped=False):
        self._filePath = filePath
        self._gzipped = gzipped

        super(EditableTerrainTile, self).fromFile(filePath, hasLighting, gzipped)

    def save(self):
        """
        persists the current state of the tile, no matter if changes were made,
        the old old state will be overwritten, if the tile is already loaded from a file
        :return: void
        """
        if not self._filePath:
            raise Exception("No _filePath defined")

        targetDirectoryPath = os.path.dirname(self._filePath)
        if not os.path.exists(targetDirectoryPath):
            os.makedirs(targetDirectoryPath)

        if os.path.exists(self._filePath):
            os.remove(self._filePath)

        self.toFile(self._filePath, self._gzipped)

    def saveTo(self, targetDirectoryPath, gzipped=False):
        """
        persists the current state of the tile into the specified directory path,
        if a tile with the same filename
        is existing, then the new file will overwrite these
        :param targetDirectoryPath: the path to the directory
        :param gzipped: whether or not the terrain tile should be gzipped
        :return: void
        """
        if self._filePath:
            tileFileName = os.path.basename(self._filePath)
        elif self.name:
            tileFileName = "{}.terrain".format(self.name)
        else:
            tileFileName = 'editable_terrain.terrain'
        if not os.path.exists(targetDirectoryPath):
            os.makedirs(targetDirectoryPath)

        filePath = os.path.join(targetDirectoryPath, tileFileName)
        if os.path.exists(filePath):
            os.remove(filePath)

        self.toFile(filePath, gzipped)

    def toWKT(self, file_path):
        """
        for debug use. persists the tile data as wkt data, all vertices and triangles
        will be create as WGS84 POINT Z and POLYGON Z WKT-Strings
        :param file_path: the file path where the wkt should be written
        :return:void
        """

        if self.isIndexDirty:
            self._rebuildIndices()

        with open(file_path, mode='w') as stream:
            vertices = self.getVerticesCoordinates()
            for i in range(len(vertices)):
                v = vertices[i]
                stream.write("POINT Z( {0} {1} {2}); {3}\n".format(v[0], v[1], v[2], i))

            indices = iter(self.indices)
            for i in range(0, len(self.indices) - 1, 3):
                vi1 = next(indices)
                vi2 = next(indices)
                vi3 = next(indices)
                llh1 = self._UVH2LLH(vi1)
                llh2 = self._UVH2LLH(vi2)
                llh3 = self._UVH2LLH(vi3)
                v1_str = "{:.14f} {:.14f} {:.14f}".format(llh1[0], llh1[1], llh1[2])
                v2_str = "{:.14f} {:.14f} {:.14f}".format(llh2[0], llh2[1], llh2[2])
                v3_str = "{:.14f} {:.14f} {:.14f}".format(llh3[0], llh3[1], llh3[2])

                stream.write("POLYGON Z(({0},{1},{2},{0})); {3}\n".format(v1_str, v2_str,
                                                                          v3_str, i))

    def findAndSplitTriangle(self, indexOfPreviousVertex, indexOfNextVertex,
                             splittingCoordinate):
        """
        Finds and splits the triangle, specified by the indexOfPreviousVertex and
        indexOfNextVertex into two new triangles with splittingCoordinate as new
        vertex of both triangles
        :param indexOfPreviousVertex:the index of the previous vertex for the new vertex
        :param indexOfNextVertex:the index of the next vertex for the new vertex
        :param splittingCoordinate: the wgs84 coordinate of the vertex between
                previous vertex and next vertex
        :return: the index of the new vertex
        """

        triangleIndex = self.findTriangleOf(indexOfPreviousVertex, indexOfNextVertex)
        if triangleIndex is None:
            raise Exception('No triangle found for Vertex')

        self.isIndexDirty = True
        oldTriangle = list(self.getTriangleAt(triangleIndex))
        newTriangle = list(oldTriangle)

        longitude, latitude, height = splittingCoordinate
        u = self._quantizeLongitude(longitude)
        v = self._quantizeLatitude(latitude)

        # insert new vertex in u,v,h
        self.u.append(u)
        self.v.append(v)
        indexOfNewVertex = len(self.u) - 1

        if self.header['minimumHeight'] < height < self.header['maximumHeight']:
            if self._changedHeights:
                self._changedHeights.append(height)
            h = self._quantizeHeight(height)
        else:
            if not self._changedHeights:
                self._changedHeights = [self._dequantizeHeight(x) for x in self.h]
            self._changedHeights.append(height)
            h = 0

        self.h.append(h)

        if type(self.vLight) == numpy.ndarray:
            self.vLight = self.vLight.tolist()

        self.vLight.append(null_normal)

        # update triangle with new vertex index
        vertexOffset = oldTriangle.index(indexOfNextVertex)
        oldTriangle[vertexOffset] = indexOfNewVertex

        # create new triangle with 'vertex_insert'
        newTriangle[newTriangle.index(indexOfPreviousVertex)] = indexOfNewVertex

        triangleOffset = (triangleIndex * 3)
        # update old triangle in indices-Array
        self.indices[int(triangleOffset + vertexOffset)] = indexOfNewVertex
        # add new triangle to indices-Array
        if type(self.indices) == numpy.ndarray:
            self.indices = list(self.indices)
        self.indices.extend(newTriangle)
        self._triangles = []

        return indexOfNewVertex

    def rebuildH(self):
        """
        Requantize the heights and sets min/max heights of this tile, if heights are
        changed, otherwise nothing will happens
        """
        if self._changedHeights:
            newMax = max(self._changedHeights)
            newMin = min(self._changedHeights)

            deniv = newMax - newMin
            b_height = old_div(TerrainTile.MAX, deniv)
            for i in range(len(self._changedHeights)):
                changed_height = self._changedHeights[i]
                h = int(round((changed_height - newMin) * b_height))
                if h < 0:
                    h = 0
                if h > TerrainTile.MAX:
                    h = TerrainTile.MAX
                self.h[i] = h

            self.header['minimumHeight'] = newMin
            self.header['maximumHeight'] = newMax
            self._changedHeights = []

    def _rebuildIndices(self):
        """
        Private method, should only be used internally if any edits
        on self.u, self.v, self.h  are made.
        """
        if self._changedHeights:
            self.rebuildH()
        sizeIndices = len(self.indices)
        sizeUVH = len(self.u)
        newU = [None] * sizeUVH
        newV = [None] * sizeUVH
        newH = [None] * sizeUVH
        newVLight = [None] * sizeUVH
        newIndices = [None] * sizeIndices
        indexMap = [None] * sizeUVH

        newIndex = 0
        for position, oldIndex in enumerate(self.indices):
            if indexMap[oldIndex]:
                (index, positions) = indexMap[oldIndex]
                positions.append(position)
            else:
                indexMap[oldIndex] = (newIndex, [position])
                newIndex += 1

        for oldIndex, data in enumerate(indexMap):
            (newIndex, positions) = data
            newU[newIndex] = (self.u[oldIndex])
            newV[newIndex] = (self.v[oldIndex])
            newH[newIndex] = (self.h[oldIndex])
            newVLight[newIndex] = (self.vLight[oldIndex])

            for position in positions:
                newIndices[position] = newIndex

        if len(self.indices) == len(newIndices):
            self.indices = newIndices
        else:
            raise Exception("Array-Size of Indices not equal")

        self.u = newU
        self.v = newV
        self.h = newH
        self.vLight = newVLight

        self.westI = self.getEdgeIndices('w')
        self.southI = self.getEdgeIndices('s')
        self.eastI = self.getEdgeIndices('e')
        self.northI = self.getEdgeIndices('n')

    def _quantizeLatitude(self, latitude):
        """
        Private helper method to convert latitude values to quantized tile (v) values
        :param latitude: the wgs 84 latitude in degrees
        :return: the quantized value (v)
        """
        b_lat = old_div(TerrainTile.MAX, (self._north - self._south))
        v = int(round((latitude - self._south) * b_lat))
        return v

    def _quantizeLongitude(self, longitude):
        """
        Private helper method to convert longitude values to quantized tile (u) values
        :param longitude: the wgs 84 longitude in degrees
        :return: the quantized value (u)
        """
        b_lon = old_div(TerrainTile.MAX, (self._east - self._west))
        u = int(round((longitude - self._west) * b_lon))
        return u

    def _quantizeHeight(self, height):
        """
        Private helper method to convert height values to quantized tile (h) values
        :param height: the wgs 84 height in ground units (meter)
        :return: the quantized value (h)
        """
        deniv = self.header['maximumHeight'] - self.header['minimumHeight']
        # In case a tile is completely flat
        if deniv == 0:
            h = 0
        else:
            b_height = old_div(TerrainTile.MAX, deniv)
            h = int(round((height - self.header['minimumHeight']) * b_height))
        return h

    def _dequantizeHeight(self, h):
        """
        Private helper method to convert quantized tile (h) values to real world height
        values
        :param h: the quantized height value
        :return: the height in ground units (meter)
        """
        return lerp(self.header['minimumHeight'],
                    self.header['maximumHeight'],
                    old_div(float(h), TerrainTile.MAX))

    def _UVH2LLH(self, index):
        """
        Private helper method to convert quantized tile vertex to wgs84 coordinate
        :param index: the index of the specified vertex
        :return: wgs84 coordinate
        """
        longitude = (
            lerp(self._west, self._east, old_div(float(self.u[index]), TerrainTile.MAX)))
        latitude = (lerp(self._south, self._north,
                         old_div(float(self.v[index]), TerrainTile.MAX)))
        height = self._dequantizeHeight(self.h[index])
        return longitude, latitude, height
