# CMakeList.txt : CMake file to build required libraries
#########
# Boost #
#########
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS regex iostreams)

# Find OpenMP
find_package(OpenMP REQUIRED)

##########
# LIBSVM #
##########
# Add Library
add_library(libsvm 
    STATIC ${CMAKE_SOURCE_DIR}/lib/svm.cpp 
    ${CMAKE_SOURCE_DIR}/include/libsvm/svm.h
)

#############
# Utilities #
#############
add_library(utils 
    STATIC ${CMAKE_SOURCE_DIR}/lib/util.cpp 
    ${CMAKE_SOURCE_DIR}/include/util.hpp
)

########################
# Off-target Penalties #
########################
add_library(otScorePenalties 
    STATIC ${CMAKE_SOURCE_DIR}/lib/otScorePenalties.cpp 
    ${CMAKE_SOURCE_DIR}/include/otScorePenalties.hpp
)

##########
# LOGGER #
##########
# Add Library
add_library(Logger 
    STATIC ${CMAKE_SOURCE_DIR}/lib/Logger.cpp
    ${CMAKE_SOURCE_DIR}/include/Logger.hpp
)

########################
# Config Parser Module #
########################
add_library(configParserModule 
    STATIC ${CMAKE_SOURCE_DIR}/lib/configParserModule.cpp 
    ${CMAKE_SOURCE_DIR}/include/configParserModule.hpp 
    ${CMAKE_SOURCE_DIR}/include/pipelineModule.hpp
)
target_link_libraries(configParserModule 
    utils 
    Boost::boost 
    Boost::regex
)

################
# Input Module #
################
add_library(inputModule 
    STATIC ${CMAKE_SOURCE_DIR}/lib/inputModule.cpp 
    ${CMAKE_SOURCE_DIR}/include/inputModule.hpp 
    ${CMAKE_SOURCE_DIR}/include/pipelineModule.hpp
)
target_link_libraries(inputModule 
    utils 
    Boost::boost 
    Boost::regex
)

#####################
# Cas9 Input Module #
#####################
add_library(cas9InputModule 
    STATIC ${CMAKE_SOURCE_DIR}/lib/cas9InputModule.cpp 
    ${CMAKE_SOURCE_DIR}/include/cas9InputModule.hpp
)
target_link_libraries(cas9InputModule 
    utils 
    Boost::regex 
    inputModule
)

####################
# Consensus Module #
####################
add_library(consensusModule 
    STATIC ${CMAKE_SOURCE_DIR}/lib/consensusModule.cpp 
    ${CMAKE_SOURCE_DIR}/include/consensusModule.hpp 
    ${CMAKE_SOURCE_DIR}/include/pipelineModule.hpp
)

###################
# CHOPCHOP Module #
###################
add_library(chopchopModule 
    STATIC ${CMAKE_SOURCE_DIR}/lib/chopchopModule.cpp 
    ${CMAKE_SOURCE_DIR}/include/chopchopModule.hpp
)
target_link_libraries(chopchopModule 
    consensusModule
)

#################
# MM10DB Module #
#################
add_library(mm10dbModule 
    STATIC ${CMAKE_SOURCE_DIR}/lib/mm10dbModule.cpp 
    ${CMAKE_SOURCE_DIR}/include/mm10dbModule.hpp
)
target_link_libraries(mm10dbModule 
    consensusModule 
    Boost::boost 
    Boost::regex
)

#######################
# sgRNAScorer2 Module #
#######################
add_library(sgrnascorer2Module 
    STATIC ${CMAKE_SOURCE_DIR}/lib/sgrnascorer2Module.cpp 
    ${CMAKE_SOURCE_DIR}/include/sgrnascorer2Module.hpp
)
target_link_libraries(sgrnascorer2Module 
    consensusModule 
    libsvm
)

######################
# Specificity Module #
######################
add_library(specificityModule 
    STATIC ${CMAKE_SOURCE_DIR}/lib/specificityModule.cpp 
    ${CMAKE_SOURCE_DIR}/include/specificityModule.hpp 
    ${CMAKE_SOURCE_DIR}/include/pipelineModule.hpp
)

##################
# Bowtie2 Module #
##################
add_library(bowtie2Module 
    STATIC ${CMAKE_SOURCE_DIR}/lib/bowtie2Module.cpp 
    ${CMAKE_SOURCE_DIR}/include/bowtie2Module.hpp
)
target_link_libraries(bowtie2Module 
    specificityModule 
    Boost::boost
)

##################################
# ISSL Off-target Scoring Module #
##################################
# Add Library
add_library(ISSLScoringModule 
    STATIC ${CMAKE_SOURCE_DIR}/lib/ISSLScoringModule.cpp 
    ${CMAKE_SOURCE_DIR}/include/ISSLScoringModule.hpp
)
# Link OpenMP if available
target_link_libraries(ISSLScoringModule 
    specificityModule 
    otScorePenalties 
    OpenMP::OpenMP_CXX 
    Boost::boost
)

#####################################################
# ISSL Off-target Scoring Module Memory Mapped File #
#####################################################
# Add Library
add_library(ISSLScoringModuleMMF 
    STATIC ${CMAKE_SOURCE_DIR}/lib/ISSLScoringModuleMMF.cpp 
    ${CMAKE_SOURCE_DIR}/include/ISSLScoringModuleMMF.hpp
)
# Link OpenMP if available
target_link_libraries(ISSLScoringModuleMMF 
    specificityModule 
    otScorePenalties 
    OpenMP::OpenMP_CXX 
    Boost::boost
    Boost::iostreams
    Boost::dynamic_linking
    Boost::diagnostic_definitions
)

#################
# Output Module #
#################
add_library(outputModule 
    STATIC ${CMAKE_SOURCE_DIR}/lib/outputModule.cpp 
    ${CMAKE_SOURCE_DIR}/include/outputModule.hpp 
    ${CMAKE_SOURCE_DIR}/include/pipelineModule.hpp
)