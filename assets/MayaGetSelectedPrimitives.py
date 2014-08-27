import maya.OpenMaya as OpenMaya
import re

def getPrimitiveRange(primitiveRange, primitiveIdx):
    indices = map(int, re.findall(r'(\b\d+\b)', primitiveRange))
    testArray = []
    if len(indices) == 0 or len(indices) > 2 :
        return

    idx1 = idx2 = indices[0]
    if len(indices) > 1:
        idx2 = indices[1]

    for idx in range(idx1, idx2 + 1):
        primitiveIdx.append(idx)
        testArray.append(idx)

    print primitiveRange
    print testArray


def getSelectedPrimitives():
    selectionList = OpenMaya.MSelectionList()
    OpenMaya.MGlobal.getActiveSelectionList(selectionList)
    selectedPrimitives = []
    selectionList.getSelectionStrings(selectedPrimitives)

    nPrimitives = len(selectedPrimitives)
    primitiveIdx = []
    for primitiveID in range(0, nPrimitives):
        getPrimitiveRange(selectedPrimitives[primitiveID], primitiveIdx)

    primitives = "primitives  " + ", ".join(map(str, primitiveIdx))

    print "# AUTO GENERATED "
    print "# selected primitives:"
    print "nPrimitives  " + str(len(primitiveIdx))
    print primitives


getSelectedPrimitives()

