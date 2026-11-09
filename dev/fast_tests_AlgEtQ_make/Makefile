# usage:
# make all          # writes one test file per .m file of the specs spec, specMod and specMtrx in split_tests/,
#                   # runs them in parallel and reports the ones that failed (default)
# make run          # runs again, in parallel, the test files already written
# make clean        # removes split_tests/
# make single       # emergency only: writes ../fast_tests_AlgEtQ.m, with all the tests, and runs it in a single
#                   # Magma session (the old behaviour, 5/6 minutes, without the report). It is not used by make all.
#
# The number of test files run at the same time is JOBS, by default the number of processors: e.g. make all JOBS=4.
# The parallel runs only need make itself (option -j); with JOBS=1 the test files are run one at a time.
# A test file fails if Magma exits with a non-zero status, if its log contains an error message, or if the run
# stops before the last line of the file. The logs are in split_tests/. If a test file fails, make exits with an error.

.PHONY : all generate run logs summary single clean

inputAlgEtQ = ../../spec
inputAlgEtQMod = ../../specMod
inputAlgEtQMtrx = ../../specMtrx
outfileAlgEtQ = ../fast_tests_AlgEtQ.m
splitdir = split_tests
JOBS ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

TESTS = $(wildcard $(splitdir)/test_*.m)
LOGS = $(TESTS:.m=.log)
# The slowest test files are started first, so that the total time is close to the time of the slowest one.
# The list is only used for the order: test files missing from it are run as well.
SLOW = MinimalGenerators IntermediateIdeals CRT ICMAbstract GraphOverOrders IdealsOfIndex Idl WkAbstract
SLOWLOGS = $(filter $(LOGS),$(SLOW:%=$(splitdir)/test_%.log))
ORDEREDLOGS = $(SLOWLOGS) $(filter-out $(SLOWLOGS),$(LOGS))

all:
	@start=$$(date +%s); \
	$(MAKE) --no-print-directory generate && $(MAKE) --no-print-directory run; \
	status=$$?; \
	echo "Wall time: $$(( $$(date +%s) - start )) s"; \
	exit $$status

generate: clean
	@echo "Writing the test files in $(splitdir)/..."
	@mkdir -p $(splitdir)
	@magma -b input:=$(inputAlgEtQ) prefix:=$(splitdir)/test attach_spec:=1 magma_make_fast_tests_AlgEtQ.m < /dev/null
	@magma -b input:=$(inputAlgEtQMod) prefix:=$(splitdir)/test_Mod attach_spec:=1 also_attach:=$(inputAlgEtQ),$(inputAlgEtQMtrx) magma_make_fast_tests_AlgEtQ.m < /dev/null
	@magma -b input:=$(inputAlgEtQMtrx) prefix:=$(splitdir)/test_Mtrx attach_spec:=1 also_attach:=$(inputAlgEtQ) magma_make_fast_tests_AlgEtQ.m < /dev/null
	@ls $(splitdir)/test_*.m > /dev/null 2>&1 || { echo "ERROR: no test file was written"; exit 1; }
	@echo "...Done: $$(ls $(splitdir)/test_*.m | wc -l) test files."

run:
	@rm -f $(splitdir)/test_*.log
	@# Attach the specs once, so that the package files are compiled before the parallel runs:
	@# otherwise all of them would write the same .sig files at the same time.
	@echo 'AttachSpec("$(inputAlgEtQ)"); AttachSpec("$(inputAlgEtQMod)"); AttachSpec("$(inputAlgEtQMtrx)"); quit;' | magma -b > /dev/null
	@echo "Running the test files, $(JOBS) at a time..."
	@$(MAKE) --no-print-directory -j$(JOBS) logs
	@$(MAKE) --no-print-directory summary

logs: $(ORDEREDLOGS)

# Magma reads /dev/null as input: otherwise, after an error, it would wait for input forever.
# The recipe always succeeds, so that all the test files are run: failures are detected by summary.
$(splitdir)/%.log: $(splitdir)/%.m
	@magma -b $< < /dev/null > $@ 2>&1; echo "### Magma exit status: $$?" >> $@
	@echo "  finished $(notdir $<)"

summary:
	@failed=0; \
	for log in $(LOGS); do \
		name=$$(basename $$log .log); \
		t=$$(grep -o ' time [0-9.]*' $$log | tail -1); \
		if grep -q '^### Magma exit status: 0$$' $$log \
			&& grep -q '^### All the tests in this file were run$$' $$log \
			&& ! grep -Eq '(Runtime|User|Internal|System) error|^In file ".*", line [0-9]+, column [0-9]+:' $$log; then \
			printf "  ok    %-40s%s\n" "$$name" "$$t"; \
		else \
			printf "  FAIL  %-40s see $(splitdir)/%s.log\n" "$$name" "$$name"; \
			failed=$$((failed+1)); \
		fi; \
	done; \
	if [ $$failed -gt 0 ]; then \
		echo "$$failed of $(words $(LOGS)) test files FAILED."; exit 1; \
	else \
		echo "All $(words $(LOGS)) test files passed."; \
	fi

# Emergency only, not used by make all.
single:
	@echo "Writing $(outfileAlgEtQ)..."
	-rm -f $(outfileAlgEtQ)
	magma -b input:=$(inputAlgEtQ) output_file:=$(outfileAlgEtQ) magma_make_fast_tests_AlgEtQ.m < /dev/null
	@echo "Running $(outfileAlgEtQ) in a single Magma session. It takes 5/6 minutes."
	magma -b $(outfileAlgEtQ) < /dev/null

clean:
	@rm -rf $(splitdir)
