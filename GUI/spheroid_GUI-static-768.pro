CONFIG      += uitools
 
INCLUDEPATH  += C:/VTK-VS8/include/vtk-5.10
INCLUDEPATH  += c:/qwt-5.2.1/src
INCLUDEPATH  += E:/ffmpeg/include

FORMS         += spheroid_GUI.ui \
                 SimpleView3DUI.ui \
                 SimpleView2DUI.ui
HEADERS       = mainwindow.h qmylabel.h params.h plot.h log.h myvtk.h misc.h \
                result_set.h graphs.h transfer.h libspheroid.h field.h \
                ImageSave.h SimpleView2DUI.h SimpleView3DUI.h \
                dialog.h myqgraphicsview.h qvideooutput.h qmycheckbox.h
RESOURCES     += icons.qrc
SOURCES       = main.cpp mainwindow.cpp params.cpp plot.cpp \
                myvtk.cpp misc.cpp lognormal.cpp graphs.cpp field.cpp \
                ImageSave.cpp SimpleView3DUI.cxx SimpleView2DUI.cxx \
                dialog.cpp myqgraphicsview.cpp  qvideooutput.cpp

# See cmake_link_command.txt for the full list of libraries that CMake links
LIBS += -LC:/VTK-VS8/lib/vtk-5.10 -lQVTK -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkIO -lvtkFiltering -lvtkCommon \
-lvtkpng -lvtktiff -lvtkjpeg -lvtkexpat -lvfw32 -lopengl32  \
-lwsock32 -lvtksys -lws2_32 -lvtkexoIIc -lvtkNetCDF \
-lvtklibxml2 -lvtkzlib -lvtkalglib \
-lgdi32 -lkernel32 -luser32 -lgdi32 -lwinspool -lshell32 -lole32 -loleaut32 -luuid -lcomdlg32 -ladvapi32    #-lpthread

LIBS += -Lc:/qwt-5.2.1/lib -lqwt5
LIBS += C:/bin/spheroid.lib

QMAKE_LIBDIR += E:/ffmpeg/lib
win32:LIBS += avcodec.lib
win32:LIBS += avdevice.lib
win32:LIBS += avfilter.lib
win32:LIBS += avformat.lib
win32:LIBS += avutil.lib
win32:LIBS += postproc.lib
win32:LIBS += swresample.lib
win32:LIBS += swscale.lib

QT           += network

#QMAKE_CXXFLAGS += -Wno-deprecated -Wno-write-strings
#NOTE: fix for link with Qt-4.7
#QMAKE_LFLAGS    += -Wl,-enable-auto-import

# install
target.path = .
sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS icons 
sources.path = .
INSTALLS += target sources

DEFINES += __MSVS_DLL__
DEFINES += _CRT_SECURE_NO_WARNINGS
DEFINES += __DISPLAY768

QMAKE_LFLAGS += /OPT:NOREF
