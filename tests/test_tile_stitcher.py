# -*- coding: utf-8 -*-
import os
import unittest

from quantized_mesh_tile import TerrainTile, tile_stitcher
from quantized_mesh_tile.tile_stitcher import TileStitcher, EdgeConnection
from tests import data_utils


class TestTileStitcher(unittest.TestCase):

    def setUp(self):
        self.quantizedTriangles = data_utils.readQuantizedTriangles()

    def testConstructor(self):
        # arrange
        centerX = 17388
        centerY = 12517
        centerZ = 14

        neighbourX = 17388
        neighbourY = 12518
        neighbourZ = 14

        # act
        centerTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                 centerX,
                                                 centerY,
                                                 centerZ)
        neighbourTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                    neighbourX,
                                                    neighbourY,
                                                    neighbourZ)
        TileStitcher(centerTile)

        # assert
        self.assertIsInstance(centerTile, TerrainTile)
        self.assertIsInstance(neighbourTile, TerrainTile)

    def testGetEdgeConnection(self):
        # arrange
        centerX = 17388
        centerY = 12517
        centerZ = 14

        neighbourX = 17388
        neighbourY = 12518
        neighbourZ = 14

        # act
        centerTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                 centerX,
                                                 centerY,
                                                 centerZ)
        neighbourTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                    neighbourX,
                                                    neighbourY,
                                                    neighbourZ)
        stitcher = TileStitcher(centerTile)
        edgeConnection = stitcher._getEdgeConnection(neighbourTile)

        # assert
        self.assertIs(edgeConnection, 'n')
        self.assertIsNotNone(edgeConnection)

    def testStitchTogetherWithSouth(self):
        # arrange
        centerX = 4347
        centerY = 3128
        centerZ = 12

        neighbourX = 4347
        neighbourY = 3127
        neighbourZ = 12

        centerTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                 centerX,
                                                 centerY,
                                                 centerZ)
        neighbourTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                    neighbourX,
                                                    neighbourY,
                                                    neighbourZ)

        # act
        stitcher = TileStitcher(centerTile)
        stitcher.addNeighbour(neighbourTile)
        stitcher.stitchTogether()
        stitcher.saveTo(data_utils.getTempPath())

        # assert
        centerTile = tile_stitcher.loadTile(
            os.path.join(data_utils.getTempPath(), '12_4347_3128.terrain'),
            centerX,
            centerY,
            centerZ)
        neighbourTile = tile_stitcher.loadTile(
            os.path.join(data_utils.getTempPath(), '12_4347_3127.terrain'),
            neighbourX,
            neighbourY,
            neighbourZ)

        centerVerticesCount = len(centerTile.getEdgeIndices(edge='s'))
        neighbourVerticesCount = len(neighbourTile.getEdgeIndices(edge='n'))

        self.assertTrue(centerVerticesCount == neighbourVerticesCount)

    def testStitchWithWestEast(self):
        # arrange
        centerX = 4347
        centerY = 3128
        centerZ = 12

        neighbourX = 4348
        neighbourY = 3128
        neighbourZ = 12

        centerTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                 centerX,
                                                 centerY,
                                                 centerZ)
        neighbourTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                    neighbourX,
                                                    neighbourY,
                                                    neighbourZ)

        # act
        stitcher = TileStitcher(centerTile)
        stitcher.addNeighbour(neighbourTile)
        stitcher.stitchTogether()

        # assert
        centerVerticesCount = len(centerTile.getEdgeIndices(edge='e'))
        neighbourVerticesCount = len(neighbourTile.getEdgeIndices(edge='w'))
        self.assertTrue(centerVerticesCount == neighbourVerticesCount)

    def testStitchWithEastAndSouth(self):
        # arrange
        centerX = 4346
        centerY = 3127
        centerZ = 12

        eastX = 4347
        eastY = 3127
        eastZ = 12

        southX = 4346
        southY = 3126
        southZ = 12

        centerTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                 centerX,
                                                 centerY,
                                                 centerZ)
        eastTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                               eastX,
                                               eastY,
                                               eastZ)
        southTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                southX,
                                                southY,
                                                southZ)

        # act
        stitcher = TileStitcher(centerTile)
        stitcher.addNeighbour(eastTile)
        stitcher.addNeighbour(southTile)
        stitcher.stitchTogether()
        stitcher.saveTo(data_utils.getTempPath())

        # assert
        centerToEastVerticesCount = len(centerTile.getEdgeIndices(edge='e'))
        centerToSouthVerticesCount = len(centerTile.getEdgeIndices(edge='s'))
        eastVerticesCount = len(eastTile.getEdgeIndices(edge='w'))
        southVerticesCount = len(southTile.getEdgeIndices(edge='n'))

        self.assertTrue(centerToEastVerticesCount == eastVerticesCount)
        self.assertTrue(centerToSouthVerticesCount == southVerticesCount)

    def testStitchWithEastAndSouth_z14x17380y12516(self):
        # arrange
        centerX = 17380
        centerY = 12516
        centerZ = 14

        eastX = 17381
        eastY = 12516
        eastZ = 14

        southX = 17380
        southY = 12515
        southZ = 14

        centerTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                 centerX,
                                                 centerY,
                                                 centerZ)
        eastTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                               eastX,
                                               eastY,
                                               eastZ)
        southTile = data_utils.buildTerrainTile(self.quantizedTriangles,
                                                southX,
                                                southY,
                                                southZ)

        # act
        stitcher = TileStitcher(centerTile)
        stitcher.addNeighbour(eastTile)
        stitcher.addNeighbour(southTile)
        stitcher.stitchTogether()

        # assert
        centerToEastVerticesCount = len(centerTile.getEdgeIndices(edge='e'))
        centerToSouthVerticesCount = len(centerTile.getEdgeIndices(edge='s'))
        eastVerticesCount = len(eastTile.getEdgeIndices(edge='w'))
        southVerticesCount = len(southTile.getEdgeIndices(edge='n'))

        self.assertTrue(centerToEastVerticesCount == eastVerticesCount)
        self.assertTrue(centerToSouthVerticesCount == southVerticesCount)

    def testHarmonizeWithEastAndSouth(self):
        # arrange
        expectedChanges = 71
        centerX = 67
        centerY = 49
        centerZ = 6

        eastX = 68
        eastY = 49
        eastZ = 6

        southX = 67
        southY = 48
        southZ = 6

        centerTile = data_utils.getTile(centerZ, centerX, centerY)
        eastTile = data_utils.getTile(eastZ, eastX, eastY)
        southTile = data_utils.getTile(southZ, southX, southY)
        normalVectorsBefore = list(centerTile.vLight)

        # act
        stitcher = TileStitcher(centerTile)
        stitcher.addNeighbour(eastTile)
        stitcher.addNeighbour(southTile)
        stitcher.harmonizeNormals()
        normalVectorsAfter = list(centerTile.vLight)

        # assert
        changes = []
        for i, normal in enumerate(normalVectorsBefore):
            normalsBefore = set(normal)
            normalsAfter = set(normalVectorsAfter[i])
            changed = normalsBefore.difference(normalsAfter)
            if changed:
                changes.append(i)

        actualChanges = len(changes)
        self.assertTrue(actualChanges == expectedChanges)

    def testGetNeighbours(self):
        # arrange
        centerX = 17380
        centerY = 12516
        centerZ = 14

        expectedWest = [14, 17379, 12516]
        expectedNorth = [14, 17380, 12517]
        expectedEast = [14, 17381, 12516]
        expectedSouth = [14, 17380, 12515]

        # act
        neighbours = tile_stitcher.getNeighbours(centerZ, centerX, centerY)

        # assert
        self.assertSequenceEqual(neighbours['west'], expectedWest)
        self.assertSequenceEqual(neighbours['north'], expectedNorth)
        self.assertSequenceEqual(neighbours['east'], expectedEast)
        self.assertSequenceEqual(neighbours['south'], expectedSouth)

    def testGetNeighboursSouthEast(self):
        # arrange
        centerX = 17380
        centerY = 12516
        centerZ = 14

        expectedEast = [14, 17381, 12516]
        expectedSouth = [14, 17380, 12515]

        # act
        neighbours = tile_stitcher.getNeighboursSouthEast(centerZ, centerX, centerY)

        # assert
        self.assertSequenceEqual(neighbours['east'], expectedEast)
        self.assertSequenceEqual(neighbours['south'], expectedSouth)

    def testEdgeConnection_repr(self):
        # arrange
        expected_repr_start = 'E:w [1] -> ({'

        # act
        ec = EdgeConnection('w', 1)
        ec.addSide('w', 2)
        ec.addSide('c', 3)
        actual_repr = ec.__repr__()

        # assert
        self.assertTrue(actual_repr.startswith(expected_repr_start))
