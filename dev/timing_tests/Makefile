# Timing tests for AlgEt. See README.md.
#
# make all [JOBS=4] [SPEC=../../spec] [MACHINE="laptop aerox"] [LABEL=...] [TIMEOUT=2400]
#     runs every test in tests/, JOBS at a time, prints a summary and appends one line per test to timings.log.
# make run TESTS="PicardGroup_1.1 OverOrders_3.*"   runs a subset (globs on the test names are allowed).
# make summary                                       prints the table of the current logs/ again.
# make compare A=<commit> B=<commit> [MACHINE=...]    side-by-side medians from timings.log (same machine).
# make list | make clean

SHELL := /bin/bash

JOBS ?= 4
SPEC ?= ../../spec
TIMEOUT ?= 2400
MACHINE ?= $(shell hostname)
LABEL ?= -
TESTS ?=
# The reference release: a timeout of this commit is never a regression (see the target log).
BASELINE ?= v1.0.2

TESTDIR := tests
LOGDIR ?= logs
LOGFILE ?= timings.log
REGRESSIONS ?= regressions.txt
# The dates in these two files are CET on every machine: UTC+1 all year, no summer time.
LOGTZ := CET-1

# The commit of the tree that SPEC points to: the exact tag or the short hash of a git checkout (plus -dirty if there
# are uncommitted changes outside this folder), or the content of the file COMMIT next to the spec of a `git archive`.
SPECDIR := $(dir $(abspath $(SPEC)))
COMMIT ?= $(shell cd $(SPECDIR) && if git rev-parse --git-dir > /dev/null 2>&1; then \
	c=$$(git describe --exact-match --tags 2>/dev/null || git rev-parse --short HEAD); \
	git diff --quiet HEAD -- . ':!dev/timing_tests' 2>/dev/null || c=$$c-dirty; echo $$c; \
	else cat COMMIT 2>/dev/null || echo unknown; fi)
MAGMAVER := $(shell echo 'a,b,c:=GetVersion(); printf "%o.%o-%o\n",a,b,c; quit;' | magma -b 2>/dev/null | tail -1)

ALLTESTS := $(sort $(wildcard $(TESTDIR)/*.m))
ifeq ($(strip $(TESTS)),)
SELECTED := $(ALLTESTS)
else
SELECTED := $(sort $(foreach t,$(TESTS),$(wildcard $(TESTDIR)/$(t).m)))
endif
LOGS := $(patsubst $(TESTDIR)/%.m,$(LOGDIR)/%.log,$(SELECTED))

# The tests that took longest on this machine (from timings.log) are started first, so that the wall time of a run
# is close to the longest test. A timeout counts as its limit (TIMEOUT(2400s) as 2400 s). Tests without a recorded
# time come last.
ORDERED := $(shell if [ -f $(LOGFILE) ]; then awk -F'\t' -v m="$(MACHINE)" '!/^\#/ && $$2==m && $$8=="ok" {s[$$6]+=$$7;n[$$6]++} \
	!/^\#/ && $$2==m && $$8 ~ /^TIMEOUT\([0-9]+s\)$$/ {x=$$8; gsub(/[^0-9]/,"",x); s[$$6]+=x;n[$$6]++} END{for(k in s) printf "%f %s\n", s[k]/n[k], k}' $(LOGFILE) | sort -rn | awk '{print "$(LOGDIR)/"$$2".log"}'; fi)
LOGS_ORDERED := $(filter $(LOGS),$(ORDERED)) $(filter-out $(ORDERED),$(LOGS))

.PHONY: all run logs summary log compare list clean check

all: clean
	@$(MAKE) --no-print-directory run

run:
	@[ -n "$(SELECTED)" ] || { echo "No test selected (TESTS='$(TESTS)')."; exit 1; }
	@mkdir -p $(LOGDIR)
	@echo "Timing tests: $(words $(SELECTED)) tests, $(JOBS) at a time, timeout $(TIMEOUT)s."
	@echo "  spec $(SPEC) (commit $(COMMIT)), machine '$(MACHINE)', Magma $(MAGMAVER), label '$(LABEL)'"
	@# Attach the spec once, so that the package is compiled before the parallel runs start.
	@echo 'AttachSpec("$(SPEC)"); quit;' | magma -b > /dev/null
	@start=$$(date +%s); $(MAKE) --no-print-directory -j$(JOBS) logs; \
	 echo "Wall time of the run: $$(( $$(date +%s) - start )) s"
	@$(MAKE) --no-print-directory summary
	@$(MAKE) --no-print-directory log

logs: $(LOGS_ORDERED)

# Magma reads /dev/null as input: otherwise, after an error, it would wait for input forever.
# The recipe always succeeds, so that all the tests are run; failures and timeouts are read from the logs.
$(LOGDIR)/%.log: $(TESTDIR)/%.m
	@mkdir -p $(LOGDIR)
	@start=$$(date +%s); \
	 timeout $(TIMEOUT) magma -b spec:=$(SPEC) name:=$* $< < /dev/null > $@ 2>&1; st=$$?; \
	 echo "### EXIT $$st" >> $@; echo "### WALL $$(( $$(date +%s) - start ))" >> $@; echo "### LIMIT $(TIMEOUT)" >> $@; \
	 t=$$(grep -m1 '^### TIMING' $@ | awk '{print $$4}'); \
	 if [ "$$st" = "124" ]; then r="TIMEOUT after $(TIMEOUT)s"; elif [ "$$st" = "0" ] && [ -n "$$t" ]; then r="$$t s"; else r="FAILED (exit $$st)"; fi; \
	 echo "  finished $* : $$r"

# status: ok (exit 0, a TIMING line and the DONE line), TIMEOUT(<limit>s) (timeout exit status 124; the limit in
# seconds is the one written in the log by the run), FAIL otherwise.
define READ_LOG
	name=$$(basename $$log .log); \
	st=$$(grep -m1 '^### EXIT' $$log | awk '{print $$3}'); \
	t=$$(grep -m1 '^### TIMING' $$log | awk '{print $$4}'); \
	w=$$(grep -m1 '^### WALL' $$log | awk '{print $$3}'); \
	lim=$$(grep -m1 '^### LIMIT' $$log | awk '{print $$3}'); lim=$${lim:-$(TIMEOUT)}; \
	if [ "$$st" = "124" ]; then status="TIMEOUT($${lim}s)"; t=-; \
	elif [ "$$st" = "0" ] && grep -q '^### DONE' $$log && [ -n "$$t" ]; then status=ok; \
	else status=FAIL; t=$${t:--}; fi
endef

summary:
	@printf "%-36s %12s %8s  %s\n" "test" "cpu seconds" "wall" "status"
	@for log in $(LOGS); do [ -f $$log ] || continue; $(READ_LOG); printf "%-36s %12s %8s  %s\n" "$$name" "$$t" "$$w" "$$status"; done

# One line per test: date | machine | magma | commit | label | test | cpu seconds | status
# A timeout is logged with status TIMEOUT(<limit>s). With a commit other than BASELINE, a timeout of a test that
# completed before on this machine (with any commit) is also reported as a regression (regressions.txt).
log:
	@[ -f $(LOGFILE) ] || printf "# date\tmachine\tmagma\tcommit\tlabel\ttest\tseconds\tstatus\n" > $(LOGFILE)
	@d=$$(TZ=$(LOGTZ) date '+%Y-%m-%d %H:%M'); for log in $(LOGS); do [ -f $$log ] || continue; $(READ_LOG); \
	 printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$$d" "$(MACHINE)" "$(MAGMAVER)" "$(COMMIT)" "$(LABEL)" "$$name" "$$t" "$$status" >> $(LOGFILE); \
	 if [ "$${status%%(*}" = "TIMEOUT" ] && [ "$(COMMIT)" != "$(BASELINE)" ]; then \
	   prev=$$(awk -F'\t' -v m="$(MACHINE)" -v t="$$name" '!/^#/ && $$2==m && $$6==t && $$8=="ok" {printf "%s %ss; ", $$4, $$7}' $(LOGFILE)); \
	   if [ -n "$$prev" ]; then \
	     echo "  REGRESSION: $$name did not finish within $${lim}s on '$(MACHINE)' with $(COMMIT), but it completed before ($$prev). Investigate."; \
	     [ -f $(REGRESSIONS) ] || printf "# date\tmachine\tcommit\ttest\tearlier completed runs (commit seconds) on the same machine\n" > $(REGRESSIONS); \
	     printf "%s\t%s\t%s\t%s\t%s\n" "$$d" "$(MACHINE)" "$(COMMIT)" "$$name" "$$prev" >> $(REGRESSIONS); \
	   fi; \
	 fi; done
	@echo "Appended $(words $(LOGS)) lines to $(LOGFILE)."

# Medians over the runs of A and B on MACHINE. A and B are an exact tag or a short hash, as written in the commit
# column (a tree with uncommitted changes outside dev/timing_tests is logged as <commit>-dirty and matches too).
# A tree whose runs on a test all timed out shows >limit, and the change is then a bound (rounded toward 0): <-65%
# means at least 65% faster, >+200% at least 200% slower.
compare:
	@[ -n "$(A)" ] && [ -n "$(B)" ] || { echo "usage: make compare A=<commit> B=<commit> [MACHINE=...]"; exit 1; }
	@awk -F'\t' -v a="$(A)" -v b="$(B)" -v m="$(MACHINE)" ' \
	  function median(str,   n,i,j,x,arr,tmp) { n=split(str,arr," "); for(i=2;i<=n;i++){x=arr[i]+0;j=i-1;while(j>=1&&arr[j]+0>x){arr[j+1]=arr[j];j--}arr[j+1]=x} \
	       if(n==0) return -1; if(n%2) return arr[(n+1)/2]+0; return (arr[n/2]+arr[n/2+1])/2 } \
	  function matches(c,x) { return c==x || c==x "-dirty" } \
	  function limit(s) { sub(/^TIMEOUT\(/,"",s); sub(/s\)$$/,"",s); return (s ~ /^[0-9]+$$/) ? s+0 : 0 } \
	  function show(x,l,h) { if (x>=0) return sprintf("%.2f",x); if (l>0) return ">" l; if (h) return "timeout"; return "-" } \
	  !/^#/ && $$2==m { if (matches($$4,a)) s="a"; else if (matches($$4,b)) s="b"; else next; names[$$6]=1; \
	    if ($$8=="ok") { if (s=="a") ta[$$6]=ta[$$6] " " $$7; else tb[$$6]=tb[$$6] " " $$7 } \
	    else if ($$8 ~ /^TIMEOUT/) { l=limit($$8); if (s=="a") { ha[$$6]=1; if (l>la[$$6]) la[$$6]=l } else { hb[$$6]=1; if (l>lb[$$6]) lb[$$6]=l } } } \
	  END { printf "%-36s %10s %10s %8s   (machine: %s)\n","test",a,b,"change",m; \
	        for (k in names) { x=median(ta[k]); y=median(tb[k]); c=""; \
	          if (x>=0 && y>=0) c=sprintf("%+.0f%%",(x>0?100*(y-x)/x:0)); \
	          else if (x<0 && la[k]>0 && y>=0) c=sprintf("<%+d%%",int(100*(y-la[k])/la[k])); \
	          else if (x>0 && y<0 && lb[k]>0) c=sprintf(">%+d%%",int(100*(lb[k]-x)/x)); \
	          printf "%-36s %10s %10s %8s\n",k,show(x,la[k],ha[k]),show(y,lb[k],hb[k]),c } }' $(LOGFILE) | (read h; echo "$$h"; sort)

list:
	@for f in $(ALLTESTS); do echo "$$(basename $$f .m)"; sed -n '2,3p' $$f | sed 's/^/      /'; done

# Checks that every test file has the required lines, does not call Random and attaches only spec.
check:
	@ok=1; for f in $(ALLTESTS); do n=$$(basename $$f .m); \
	  grep -q 'AttachSpec(spec)' $$f || { echo "$$n: missing AttachSpec(spec)"; ok=0; }; \
	  grep -q '### TIMING' $$f || { echo "$$n: missing the TIMING line"; ok=0; }; \
	  grep -q '### DONE' $$f || { echo "$$n: missing the DONE line"; ok=0; }; \
	  grep -q '^quit;' $$f || { echo "$$n: missing quit;"; ok=0; }; \
	  ! grep -q 'Random(' $$f || { echo "$$n: calls Random (the inputs must be fixed, see README)"; ok=0; }; \
	  ! grep -q 'AttachSpec(spec cat' $$f || { echo "$$n: attaches a package other than spec (see README)"; ok=0; }; \
	  grep -q "^ *$$n" $$f || { echo "$$n: the header does not start with the test name"; ok=0; }; done; \
	 [ $$ok = 1 ] && echo "All $(words $(ALLTESTS)) test files look fine."

clean:
	@rm -rf $(LOGDIR)
