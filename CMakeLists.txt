cmake_minimum_required (VERSION 3.4)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_HOST_SYSTEM "Windows-10.0.22621")
set(CMAKE_CONFIGURATION_TYPES Debug Release)

project(ArduinoLEDTester)
math(EXPR BITS "8*${CMAKE_SIZEOF_VOID_P}")

FUNCTION(ADD_SOURCE_DIRS SourcesVar)
	foreach(source IN LISTS ${SourcesVar})
		get_filename_component(source_path "${source}" PATH)
		
		if (source_path STREQUAL "")
		else (source_path STREQUAL "")
		string(FIND ${source_path} "include" foundInclude REVERSE)
		string(FIND ${source_path} "src" foundSrc REVERSE)
		string(FIND ${source_path} "thirdparty" foundThirdparty)
		if (foundInclude GREATER 0)
			string(SUBSTRING ${source_path} ${foundInclude} -1 source_path)
		endif ()
		if (foundSrc GREATER 0)
			string(SUBSTRING ${source_path} ${foundSrc} -1 source_path)
		endif ()
		if (foundThirdparty GREATER 0)
			string(SUBSTRING ${source_path} ${foundThirdparty} -1 source_path)
		endif ()
		
		string(REPLACE "/" "\\" source_path_msvc "${source_path}")
		string(REPLACE "include" "Header Files" source_path_msvc "${source_path_msvc}")
		string(REPLACE "src" "Source Files" source_path_msvc "${source_path_msvc}")
		source_group("${source_path_msvc}" FILES "${source}")
		endif (source_path STREQUAL "")
	endforeach()
ENDFUNCTION()


add_definitions(-DPLATFORM_WINDOWS)
add_definitions(-D_HAS_EXCEPTIONS=0 -D_STATIC_CPPLIB)
add_definitions(-DUNICODE -D_UNICODE)
add_definitions(-DPSAPI_VERSION=1) # https://msdn.microsoft.com/en-us/library/windows/desktop/ms682623(v=vs.85).aspx
add_definitions(-D_SILENCE_CXX17_CODECVT_HEADER_DEPRECATION_WARNING)


set(common_build_flags "/MP32 /arch:SSE2 /std:c++20 /std:c11 /Zc:tlsGuards-")

set(CMAKE_CXX_FLAGS_RELEASE "/MT /GF /GR- /Gy /Oi /O2 /Ot /Zi /GL /GS- ${common_build_flags} ") # /fp:fast

set(CMAKE_CXX_FLAGS_DEBUG "/MTd /GR- /Od /Zi /RTC1 ${common_build_flags} ")
set(CMAKE_C_FLAGS_RELEASE "/MT /GF /GR- /Gy /Oi /O2 /Ot /LTCG /Zi /GL /GS- ${common_build_flags} ")
set(CMAKE_C_FLAGS_DEBUG "/MTd /GR- /Od /Zi /RTC1 ${common_build_flags} ")

set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /INCREMENTAL:NO /DEBUG /LTCG ")
set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS_RELEASE} ")
set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /INCREMENTAL:NO /DEBUG /LTCG ")
string(REPLACE "/EHsc" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REPLACE "/EHsc" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

set(app_src
	"main.cpp"
	
	"app.manifest"
)

include_directories(".")
include_directories("include")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE /SUBSYSTEM:WINDOWS ${LINKER_DELAY_LOADING} ")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /OPT:REF /OPT:ICF /NODEFAULTLIB:\"LIBCMT\"")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /SAFESEH:NO /NODEFAULTLIB:\"LIBCMT\"")

set(projoutputname "app")
add_definitions(-DCMAKE_BINARY_NAME="${projoutputname}")
add_executable(${projoutputname} ${app_src})
target_link_options(${projoutputname} PRIVATE /PDBALTPATH:$<TARGET_PDB_FILE_NAME:${projoutputname}>)
target_link_libraries(${projoutputname} optimized libcmt)

ADD_SOURCE_DIRS(app_src)
