# usage:
# make all	    	# cleans up old file and makes a new one

.PHONY : clean_old 	# prevent errors from rm

inputAlgEtQ = ../../spec
outfileAlgEtQ = ../fast_tests_AlgEtQ.m

all: clean_old magma_script

clean_old:
	@echo "Removing old fast_tests_AlgEtQ.m..."
	-rm ${outfileAlgEtQ}
	@echo "...Done"

magma_script:
	@echo "Printing the new fast_tests_AlgEtQ.m"
	@echo "The tests will be sorted in the same way as the files they are in are listed in the spec file"
	magma -b input:=${inputAlgEtQ} output_file:=${outfileAlgEtQ} magma_make_fast_tests_AlgEtQ.txt
	@echo "...Done"
	@echo "--------------------------------------"
	@echo "--------------------------------------"
	@echo "We now run the fast tests. It should take only a few minutes."
	@echo "--------------------------------------"
	@echo "--------------------------------------"
	magma -b ../fast_tests_AlgEtQ.m