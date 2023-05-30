# usage:
# make all	    	# cleans up old file and makes a new one


SOURCES = $(pathsubst %, ../../%, $(shell grep \\.m$$ ../../spec))


inputAlgEtQ = ../../spec
outfileAlgEtQ = ../fast_tests_AlgEtQ.m

all: magma_script

clean:
	@echo "Removing old fast_tests_AlgEtQ.m..."
	if [ -f $(outfileAlgEtQ) ]; then rm $(outfileAlgEtQ) ; fi
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

print-%:
	@echo '$*=$($*)'

.PHONY : clean_old 	# prevent errors from rm

