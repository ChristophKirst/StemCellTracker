
vImarisApplication = imarisinstance;

aSurface = vImarisApplication.GetFactory.ToSurfaces(vImarisApplication.GetSurpassSelection);
vSurfaceListValues = aSurface.GetSurfacesList(0);
vVertices = vSurfaceListValues.mVertices;
vNumberOfVerticesPerSurface = vSurfaceListValues.mNumberOfVerticesPerSurface;
vTriangles = vSurfaceListValues.mTriangles;
vNumberOfTrianglesPerSurface = vSurfaceListValues.mNumberOfTrianglesPerSurface;
vNormals = vSurfaceListValues.mNormals;
vTimeIndexPerSurface = vSurfaceListValues.mTimeIndexPerSurface;