cmake_minimum_required(VERSION 3.5)

project(dimer_generator VERSION 0.1 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_COMPILE_PDB_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/pdb)
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})

# Setup for work with RDKit C++ API
set(RDBASE $ENV{RDBASE})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${RDBASE} ${RDBASE}/Code/cmake/Modules)
add_definitions("-DRDK_BUILD_CORDGEN_SUPPORT=ON")
add_compile_options(-Wall)

# set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
# set(Boost_USE_STATIC_RUNTIME ON)

find_package( QT NAMES Qt6 Qt5 REQUIRED
    COMPONENTS
        Widgets
)
find_package( Qt${QT_VERSION_MAJOR} REQUIRED
    COMPONENTS
        Widgets
        OpenGL
)
find_package( Boost REQUIRED
    COMPONENTS
        iostreams
        filesystem
        system
        serialization
        regex
)
find_package( Cairo REQUIRED )
find_package( Threads REQUIRED )
find_package( rdkit CONFIG REQUIRED )
find_package( Freetype REQUIRED )

set(RDKit_LIBS
    RDKit::FileParsers
    RDKit::SmilesParse
    RDKit::Depictor
    RDKit::RDGeneral
    RDKit::MolDraw2D
    RDKit::MolDraw2DQt
    RDKit::GraphMol
    RDKit::ChemReactions
    RDKit::MolStandardize
)

set( LIBS
    ${Boost_IOSTREAMS_LIBRARY}
    ${Boost_REGEX_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_SERIALIZATION_LIBRARY}
    ${CAIRO_LIBRARIES} )

set(PROJECT_SOURCES
    main.cpp
    mainwindow.cpp
    reaction.cpp
    molecule.cpp
    chemicalitem.cpp
    reactiondialog.cpp
    reactiongenerator.cpp
)

set(PROJECT_HEADERS
    mainwindow.h
    reaction.h
    molecule.h
    chemicalitem.h
    reactiondialog.h
    reactiongenerator.h
)

set(PROJECT_FORMS
    mainwindow.ui
    chemicalitem.ui
    reactiondialog.ui
)

set(PROJECT_RESOURCES
    res/QtResources.qrc
)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(dimer_generator
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
        ${PROJECT_HEADERS}
        ${PROJECT_FORMS}
        ${PROJECT_RESOURCES}
    )
# Define target properties for Android with Qt 6 as:
#    set_property(TARGET dimer_generator APPEND PROPERTY QT_ANDROID_PACKAGE_SOURCE_DIR
#                 ${CMAKE_CURRENT_SOURCE_DIR}/android)
# For more information, see https://doc.qt.io/qt-6/qt-add-executable.html#target-creation
else()
    if(ANDROID)
        add_library(dimer_generator SHARED
            ${PROJECT_SOURCES}
            ${PROJECT_HEADERS}
            ${PROJECT_FORMS}
            ${PROJECT_RESOURCES}
        )
# Define properties for Android with Qt 5 after find_package() calls as:
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
    else()
        add_executable(dimer_generator
            ${PROJECT_SOURCES}
            ${PROJECT_HEADERS}
            ${PROJECT_FORMS}
            ${PROJECT_RESOURCES}
        )
    endif()
endif()

target_include_directories(dimer_generator PUBLIC ${RDKit_INCLUDE_DIRS} ${CAIRO_INCLUDE_DIRS})

target_link_libraries( dimer_generator
    PUBLIC ${LIBS} ${RDKit_LIBS}
    PRIVATE Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::OpenGL Freetype::Freetype
)

set_target_properties(dimer_generator PROPERTIES
    MACOSX_BUNDLE_GUI_IDENTIFIER my.example.com
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(dimer_generator)
endif()

install()
