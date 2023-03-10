# usage:
# make all
# Submission of AlgEtQ 20230310

.PHONY : clean_old 	# prevent errors from rm

date != date +%Y%m%d
zip = ${date}_AlgEtQ_marseglia.zip

all: zip

zip:
	mkdir AlgEt
	mkdir AlgEt/AlgEtQ
	mkdir AlgEt/dev/
	mkdir AlgEt/dev/fast_tests_AlgEtQ_make
	mkdir AlgEt/doc
	mkdir AlgEt/examples
	cp ../../spec AlgEt/
	cp ../../AlgEtQ/*.m AlgEt/AlgEtQ/
	cp ../fast_tests_AlgEtQ_make/* AlgEt/dev/fast_tests_AlgEtQ_make/
	cp ../slow_tests_AlgEtQ.m AlgEt/dev/slow_tests_AlgEtQ.m
	cp ../../examples/ideal_class_monoid.txt AlgEt/examples/ideal_class_monoid.txt
	cp ../../doc/Makefile AlgEt/doc/
	cp ../../doc/magma_make* AlgEt/doc/
	find AlgEt/AlgEtQ/ -name *.m -exec sed -i 's/marse004\//marse004\/\n\/\/ last modified on ${date}/g' {} \;
	zip -r ${zip} AlgEt
	zip ${zip} README.md
	rm -r AlgEt

