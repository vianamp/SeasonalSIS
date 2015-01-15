CC = g++

DFLAG = -DDEBUG

CFLAGS = -O1 -O2 -O3 -Os

VTK_INCLUDE_DIR = /usr/local/include/vtk-6.2/

VTK_LIBS_DIR = /usr/local/lib/

VTK_DYN_LIBS = -lvtkImagingCore-6.2 -lvtkFiltersSources-6.2 -lvtkFiltersCore-6.2 -lvtkImagingGeneral-6.2 -lvtkCommonDataModel-6.2 -lvtkCommonCore-6.2 -lvtkCommonTransforms-6.2 -lvtkIOLegacy-6.2 -lvtkCommonExecutionModel-6.2 -lvtkIOCore-6.2 -lvtkIOImage-6.2 -lvtkCommonSystem-6.2 -lvtksys-6.2 -lvtkImagingStencil-6.2 -lvtkCommonMisc-6.2 -lvtkCommonTransforms-6.2 -lvtkCommonMath-6.2 -lvtktiff-6.2 -lvtkjpeg-6.2 -lvtkpng-6.2 -lvtkzlib-6.2 -lvtkRenderingFreeType-6.2 -lvtkFiltersGeneral-6.2

VTK_STA_LIBS = $(VTK_LIBS_DIR)libvtkImagingCore-6.2.a $(VTK_LIBS_DIR)libvtkFiltersSources-6.2.a $(VTK_LIBS_DIR)libvtkFiltersCore-6.2.a $(VTK_LIBS_DIR)libvtkImagingGeneral-6.2.a $(VTK_LIBS_DIR)libvtkCommonDataModel-6.2.a $(VTK_LIBS_DIR)libvtkCommonTransforms-6.2.a $(VTK_LIBS_DIR)libvtkCommonCore-6.2.a $(VTK_LIBS_DIR)libvtkIOLegacy-6.2.a $(VTK_LIBS_DIR)libvtkCommonExecutionModel-6.2.a $(VTK_LIBS_DIR)libvtkIOCore-6.2.a $(VTK_LIBS_DIR)libvtkIOImage-6.2.a $(VTK_LIBS_DIR)libvtkCommonSystem-6.2.a $(VTK_LIBS_DIR)libvtksys-6.2.a $(VTK_LIBS_DIR)libvtkImagingStencil-6.2.a $(VTK_LIBS_DIR)libvtkCommonMisc-6.2.a $(VTK_LIBS_DIR)libvtkCommonTransforms-6.2.a $(VTK_LIBS_DIR)libvtkCommonMath-6.2.a $(VTK_LIBS_DIR)libvtktiff-6.2.a $(VTK_LIBS_DIR)libvtkpng-6.2.a $(VTK_LIBS_DIR)libvtkjpeg-6.2.a $(VTK_LIBS_DIR)libvtkzlib-6.2.a $(VTK_LIBS_DIR)libvtkRenderingFreeType-6.2.a $(VTK_LIBS_DIR)libvtkFiltersGeneral-6.2.a

IGRAPH_INCLUDE_DIR = /usr/local/include/igraph/

IGRAPH_LIBS_DIR = /usr/local/lib/

IGRAPH_DYN_LIBS = -ligraph

IGRAPH_STA_LIBS = $(IGRAPH_LIBS_DIR)libigraph.a

debug: SeasonalSIS.cpp
	@$(CC) $(DFLAG) $(CFLAGS) SeasonalSIS.cpp -o SeasonalSIS $(VTK_DYN_LIBS) $(IGRAPH_DYN_LIBS) -I$(VTK_INCLUDE_DIR) -I$(IGRAPH_INCLUDE_DIR)

dynamic: SeasonalSIS.cpp
	@$(CC) $(CFLAGS) SeasonalSIS.cpp -o SeasonalSIS $(VTK_DYN_LIBS) $(IGRAPH_DYN_LIBS) -I$(VTK_INCLUDE_DIR) -I$(IGRAPH_INCLUDE_DIR)

static:
	@$(CC) $(CFLAGS) SeasonalSIS.cpp -o SeasonalSIS $(VTK_STA_LIBS) $(IGRAPH_STA_LIBS) -I$(VTK_INCLUDE_DIR) -I$(IGRAPH_INCLUDE_DIR)