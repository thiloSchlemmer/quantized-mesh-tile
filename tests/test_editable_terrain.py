# -*- coding: utf-8 -*-
import os
import unittest
from tests import data_utils


class TestEditableTerrainTile(unittest.TestCase):
    def setUp(self):
        self.quantizedTriangles = data_utils.readQuantizedTriangles()

    def testGetEdgeCoordinates(self):
        # arrange
        x = 17388
        y = 12517
        z = 14
        edge = 'w'

        # act
        tile = data_utils.buildTerrainTile(self.quantizedTriangles, x, y, z)
        coordinates = tile.getEdgeCoordinates(edge)

        # assert
        self.assertTrue(len(coordinates) == 2)
        self.assertTrue(len(coordinates[0]) == 3)

    def testToWKT(self):
        # arrange
        x = 17388
        y = 12517
        z = 14
        wktPath = os.path.join(data_utils.getTempPath(), 'test.wkt')

        try:
            # act
            tile = data_utils.buildTerrainTile(self.quantizedTriangles, x, y, z)
            tile.toWKT(wktPath)
            # assert
            with open(wktPath, mode='r') as wktFile:
                lines = wktFile.readlines()
            self.assertGreater(len(lines), 0)
        finally:
            os.remove(wktPath)

    def testSetHeight_withHeightUnderMin(self):
        # arrange
        x = 17388
        y = 12517
        z = 14
        expectedHeight = 0

        # act
        tile = data_utils.buildTerrainTile(self.quantizedTriangles, x, y, z)
        tile.setHeight(0, -1.0)
        tile.rebuildH()
        actualHeight = tile.h[0]

        # assert
        self.assertEquals(actualHeight, expectedHeight)

    def testSetHeight_withHeightOverMax(self):
        # arrange
        x = 17388
        y = 12517
        z = 14
        expectedHeight = 9999

        # act
        tile = data_utils.buildTerrainTile(self.quantizedTriangles, x, y, z)
        tile.setHeight(0, 9999)
        tile.rebuildH()
        actualHeight = tile.getHeightAt(0)

        # assert
        self.assertEquals(actualHeight, expectedHeight)

    def testSave_expectedException(self):
        # arrange
        x = 17388
        y = 12517
        z = 14

        # act
        tile = data_utils.buildTerrainTile(self.quantizedTriangles, x, y, z)

        # assert
        self.assertRaises(Exception, tile.save)

    def testSave_expectedExistingFile(self):
        # arrange
        x = 17388
        y = 12517
        z = 14

        try:
            # act
            tile = data_utils.buildTerrainTile(self.quantizedTriangles, x, y, z)
            testPath = os.path.join(data_utils.getTempPath(), "test.terrain")
            tile._filePath = testPath
            tile.save()

            # assert
            self.assertTrue(os.path.exists(testPath))
        finally:
            if os.path.exists(testPath):
                os.remove(testPath)
