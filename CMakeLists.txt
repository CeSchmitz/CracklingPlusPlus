# CMakeList.txt : Main CMake file for CracklingPlusPlus
cmake_minimum_required (VERSION 3.8)

# Set CXX Standard for project
set(CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Update to use target compile flags. temp solution
#if(NOT CMAKE_BUILD_TYPE)
	#set(CMAKE_BUILD_TYPE Release)
#endif()

#set(CMAKE_CXX_FLAGS "")
#set(CMAKE_CXX_FLAGS_DEBUG "-G")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Project name, version and description
project ("CracklingPlusPlus" LANGUAGES CXX C VERSION 1.0.0 DESCRIPTION "C++ version of Crackling" )

# Add source to this project's executable.
add_executable (${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/CracklingPlusPlus.cpp)

# Add library subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/lib)

find_package(Boost REQUIRED COMPONENTS regex iostreams)

# Link Libraries to main executable
target_link_libraries (
	${PROJECT_NAME} 
	Logger 
	Boost::dynamic_linking
	Boost::boost
	Boost::regex
	Boost::iostreams
	configParserModule
	cas9InputModule
	chopchopModule
	mm10dbModule
	sgrnascorer2Module
	bowtie2Module
	ISSLScoringModule
	ISSLScoringModuleMMF
	outputModule
)

# Added ISSLCreateIndex subdirectory
add_subdirectory (${CMAKE_SOURCE_DIR}/ISSLCreateIndex)

# Add tests subdirectory
# add_subdirectory (${CMAKE_SOURCE_DIR}/tests)

# TODO: Sort out testing data in single location
# Copy testing data
file(COPY ${CMAKE_SOURCE_DIR}/tests/data/ DESTINATION data/)
