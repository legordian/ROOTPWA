
find_package(Git QUIET)
if(GIT_FOUND)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} log --pretty="%H" -n1 -- ${CMAKE_SOURCE_DIR_}
		OUTPUT_VARIABLE GIT_HASH
		RESULT_VARIABLE _GIT_LOG_RETURN
		OUTPUT_STRIP_TRAILING_WHITESPACE
		)
endif()


execute_process(COMMAND hostname
	OUTPUT_VARIABLE HOSTNAME
	RESULT_VARIABLE _HOSTNAME_RETURN
	OUTPUT_STRIP_TRAILING_WHITESPACE
	)
if(_HOSTNAME_RETURN)
	unset(HOSTNAME)
endif()
unset(_HOSTNAME_RETURN)


set(USER $ENV{USER})


configure_file(${SRC} ${DST} @ONLY)
