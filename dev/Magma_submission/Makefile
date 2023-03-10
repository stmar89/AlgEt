# usage:
# make all
# Submission of AlgEtQ 20230310

.PHONY : clean_old 	# prevent errors from rm

date != date +%Y%m%d
zip = ${date}_AlgEtQ_marseglia.zip

all: zip

zip:
	@echo Adding spec file
	zip ${zip} ../../spec
	@echo Adding .m files from AlgEtQ
	zip ${zip} ../../AlgEtQ/*.m




