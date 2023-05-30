# usage:
# make all	    	# cleans up old file and makes a new one


SOURCES = $(pathsubst %, ../../%, $(shell grep \\.m$$ ../../spec))


inputAlgEtQ = ../../spec
outfileAlgEtQ = ../fast_tests_AlgEtQ.m

all: magma_script

clean:
	@echo "Removing old fast_tests_AlgEtQ.m..."
	if [ -f $(outfileAlgEtQ) ]; then rm $(outfileAlgEtQ) ; fi
	if [ -f set_quit_on_error.m ]; then rm set_quit_on_error.m ; fi
	@echo "...Done"

$(outfileAlgEtQ): $(SOURCES) magma_make_fast_tests_AlgEtQ.txt
	@echo "Printing the new fast_tests_AlgEtQ.m"
	@echo "The tests will be sorted in the same way as the files they are in are listed in the spec file"
	magma -b input:=$(inputAlgEtQ) output_file:=$(outfileAlgEtQ) $<
	@echo "...Done"



magma_script: $(outfileAlgEtQ)
	@echo "-------------------------------------------------------------"
	@echo "-------------------------------------------------------------"
	@echo "We now run the fast tests. It should take only a 5/6 minutes."
	@echo "-------------------------------------------------------------"
	@echo "-------------------------------------------------------------"
	magma -b $<

set_quit_on_error.m:
	-echo 'SetQuitOnError(true);\n' > set_quit_on_error.m

magma_script_fail_on_error: set_quit_on_error.m $(outfileAlgEtQ)
	@echo "-------------------------------------------------------------"
	@echo "-------------------------------------------------------------"
	@echo "We now run the fast tests. It should take only a 5/6 minutes."
	@echo "-------------------------------------------------------------"
	@echo "-------------------------------------------------------------"
	magma -b $^


print-%:
	@echo '$*=$($*)'

.PHONY : clean_old 	# prevent errors from rm

