# usage:
# make all	    	# cleans up old file and makes a new one

.PHONY : clean_old 	# prevent errors from rm

inputAlgEtQ = ../spec
outfileAlgEtQ = ./ListOfCommandsAlgEtQ.md

inputAlgEtQMod = ../specMod
outfileAlgEtQMod = ./ListOfCommandsAlgEtQMod.md

inputAlgEtQMtrx = ../specMtrx
outfileAlgEtQMtrx = ./ListOfCommandsAlgEtQMtrx.md

all: clean_old magma_script

clean_old:
	@echo "Removing old ListOfCommandsAlgEtQ.md..."
	-rm ${outfileAlgEtQ}
	@echo "...Done"
	@echo "Removing old ListOfCommandsAlgEtQMod.md..."
	-rm ${outfileAlgEtQMod}
	@echo "...Done"
	@echo "Removing old ListOfCommandsAlgEtQMtrx.md..."
	-rm ${outfileAlgEtQMtrx}
	@echo "...Done"

magma_script:
	@echo "Printing the new ListOfCommandsAlgEtQ.md."
	@echo "The commands will be sorted in the same way as they are listed in the spec file"
	magma -b input:=${inputAlgEtQ} output_file:=${outfileAlgEtQ} magma_make_list_of_func.txt
	@echo "...Done"
	@echo "Printing the new ListOfCommandsAlgEtQMod.md."
	@echo "The commands will be sorted in the same way as they are listed in the spec file"
	magma -b input:=${inputAlgEtQMod} output_file:=${outfileAlgEtQMod} magma_make_list_of_func.txt
	@echo "...Done"
	@echo "Printing the new ListOfCommandsAlgEtQMtrx.md."
	@echo "The commands will be sorted in the same way as they are listed in the spec file"
	magma -b input:=${inputAlgEtQMtrx} output_file:=${outfileAlgEtQMtrx} magma_make_list_of_func.txt
	@echo "...Done"
	@echo "...Don't forget to push the new changes!"
