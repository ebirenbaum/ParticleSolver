QT += core gui opengl

TARGET = particles_cuda
TEMPLATE = app

# project build directories
DESTDIR     = $$system(pwd) # target dir
BUILDDIR    = $$DESTDIR/build

MOC_DIR     = $$BUILDDIR # moc_...
RCC_DIR     = $$BUILDDIR # qrc_resources.cpp
UI_DIR      = $$BUILDDIR # ui_mainwindow.cpp

OBJECTS_DIR = $$BUILDDIR/bin # .o files

unix:!macx {
    NON_CUDA_LIBS = -lGLU
    LIBS += $$NON_CUDA_LIBS
    QMAKE_CXXFLAGS += -std=c++11
    DEFINES += LINUX    #define LINUX
}
macx {
    NON_CUDA_LIBS += -stdlib=libc++
    LIBS += $$NON_CUDA_LIBS

    QMAKE_CXXFLAGS += -stdlib=libc++
    QMAKE_CXXFLAGS += -std=c++11
    QMAKE_CXXFLAGS += -mmacosx-version-min=10.7
    QMAKE_LFLAGS += -mmacosx-version-min=10.7
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.7

    MAC_SDK  = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/
    if( exists( $$MAC_SDK) ) {
        QMAKE_MAC_SDK = macosx10.9
    }

    QMAKE_CXXFLAGS += -Wno-c++11-extensions
}


INCLUDEPATH += glm src src/ui src/constraints src/rendering src/cuda
DEPENDPATH += glm src src/ui src/constraints src/rendering src/cuda

SOURCES += src/main.cpp \
    src/ui/mainwindow.cpp \
    src/ui/view.cpp \
    src/constraints/distanceconstraint.cpp \
    src/rendering/renderer.cpp \
    src/particleapp.cpp \
    src/particlesystem.cpp \
    src/rendering/camera.cpp \
    src/rendering/orbitingcamera.cpp \
    src/rendering/actioncamera.cpp

HEADERS += src/ui/mainwindow.h \
    src/ui/view.h \
    src/constraints/constraint.h \
    src/constraints/distanceconstraint.h \
    src/rendering/renderer.h \
    src/particleapp.h \
    src/particlesystem.h \
    src/rendering/camera.h \
    src/rendering/orbitingcamera.h \
    src/debugprinting.h \
    src/rendering/actioncamera.h

FORMS += src/ui/mainwindow.ui

OTHER_FILES += res/shader.vert res/shader.frag

#TODO (Windows): If you are setting up local development on Windows (NOT Mac), comment out the following lines
win32:CONFIG(release, debug|release): LIBS += -L/course/cs123/lib/glew/glew-1.10.0/lib/release/ -lGLEW
else:win32:CONFIG(debug, debug|release): LIBS += -L/course/cs123/lib/glew/glew-1.10.0/lib/debug/ -lGLEW
else:unix: LIBS += -L/usr/local/Cellar/glew/1.11.0/lib -lGLEW

#TODO (Windows or Mac): If you are setting up local development on Windows OR Mac, fill in the correct path to your glew and uncomment the following lines:
INCLUDEPATH+=/usr/local/Cellar/glew/1.11.0/include
DEPENDPATH+=/usr/local/Cellar/glew/1.11.0/include

RESOURCES += \
    resources.qrc



################################ CUDA #################################


unix:!macx {
    # Path to cuda stuff
    CUDA_DIR = /contrib/projects/cuda5-toolkit
    CUDA_LIB = $$CUDA_DIR/lib64

    CUDA_DEF += -DCUDA_5

    # GPU architecture
    CUDA_ARCH     = sm_21 # should be able to detect this somehow instead of hardcoding

    SED_STUFF = 2>&1 | sed -r \"s/\\(([0-9]+)\\)/:\\1/g\" 1>&2
}
macx {
    # Path to cuda stuff
    CUDA_DIR = /usr/local/cuda
    CUDA_LIB = $$CUDA_DIR/lib

    CUDA_DEF += -DCUDA_7

    # GPU architecture
    CUDA_ARCH     = sm_30 # should be able to detect this somehow instead of hardcoding

    SED_STUFF = 2>&1 | sed -E \"s/\\(([0-9]+)\\)/:\\1/g\" 1>&2

    NVCCFLAGS += --std=c++11
}

if ( exists( $$CUDA_DIR/ ) ) {
    message( "Configuring for cuda..." );
    DEFINES += CUDA

# Cuda sources
CUDA_SOURCES += src/cuda/shared_variables.cu \
                src/cuda/integration.cu \
                src/cuda/solver.cu \
                src/cuda/util.cu

OTHER_FILES +=  src/cuda/wrappers.cuh \
                src/cuda/integration_kernel.cuh \
                src/cuda/solver_kernel.cuh \
                src/cuda/kernel.cuh \
                src/cuda/util.cuh \
                src/cuda/util.cu \
                src/cuda/integration.cu \
                src/cuda/solver.cu \
                src/cuda/shared_variables.cu \
                src/cuda/shared_variables.cuh \
                src/cuda/helper_cuda.h

    # Pather to header and lib files
    INCLUDEPATH += $$CUDA_DIR/include
    QMAKE_LIBDIR += $$CUDA_LIB

    # prevents warnings from code we didn't write
    QMAKE_CXXFLAGS += -isystem $$CUDA_DIR/include

    LIBS += -lcudart -lcublas -lcusparse -lcurand
    QMAKE_LFLAGS += -Wl,-rpath,$$CUDA_LIB
    NVCCFLAGS = -Xlinker -rpath,$$CUDA_LIB

    # libs used in the code
    CUDA_LIBS = $$LIBS
    CUDA_LIBS -= $$NON_CUDA_LIBS

    # Here are some NVCC flags I've always used by default.
    NVCCFLAGS     += --compiler-options -fno-strict-aliasing -use_fast_math --ptxas-options=-v

    # Prepare the extra compiler configuration (taken from the nvidia forum)
    CUDA_INC = $$join(INCLUDEPATH,' -I','-I',' ')

    cuda.commands = $$CUDA_DIR/bin/nvcc -m64 -O3 -arch=$$CUDA_ARCH -c $$NVCCFLAGS \
                    $$CUDA_INC $$CUDA_LIBS  ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT} \
                    $$CUDA_DEF $$SED_STUFF
    # nvcc error printout format ever so slightly different from gcc
    # http://forums.nvidia.com/index.php?showtopic=171651

    cuda.dependency_type = TYPE_C # there was a typo here. Thanks workmate!
    cuda.depend_command = $$CUDA_DIR/bin/nvcc -O3 -M $$CUDA_INC $$NVCCFLAGS   ${QMAKE_FILE_NAME}

    cuda.input = CUDA_SOURCES
    cuda.output = ${OBJECTS_DIR}${QMAKE_FILE_BASE}_cuda.o

    # Tell Qt that we want add more stuff to the Makefile
    QMAKE_EXTRA_COMPILERS += cuda
}
