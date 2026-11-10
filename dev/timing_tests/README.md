# Timing tests

A suite of timing tests for the package AlgEtQ (the files listed in `spec`), one test per file, run in parallel on a
limited number of cores, with the timings logged per machine, Magma version and commit. It replaces the file of slow
tests, now kept in `dev/old/slow_tests_AlgEtQ.m`. AlgEtNF and the packages of `specMod` (AlgEtQMod) and `specMtrx`
(AlgEtQMtrx) are not covered.

## Running

From this folder:

    make all                                   # all the tests, 4 at a time, against ../../spec
    make all JOBS=2 MACHINE="laptop aerox"     # fewer parallel jobs, a name for the machine in the log
    make all SPEC=/path/to/other/tree/spec     # the same tests against another checkout (e.g. a `git archive`
                                               # of a release: put the commit hash in a file COMMIT next to its spec)
    make run TESTS="PicardGroup_1.1 WKICM_4.*" # a subset
    make summary                               # the table of the last run again
    make compare A=v1.0.2 B=feature MACHINE="laptop aerox"   # medians from timings.log, side by side
    make list | make check | make clean

Variables: `JOBS` (parallel tests, default 4: keep it well below the number of cores, so that the timings do not
interfere), `SPEC` (default `../../spec`), `MACHINE` (free text, default `hostname`: use "laptop aerox",
"babbage" and "laplace OvGU"), `LABEL` (a free comment stored in the log), `TIMEOUT` (seconds per test, default 2400),
`TESTS`, `BASELINE` (the reference release, default v1.0.2, see below).

Each run appends one line per test to `timings.log` (tab separated: date, machine, Magma version, commit, label,
test, CPU seconds, status). The date is in CET (UTC+1 all year, no summer time) on every machine. The commit is the
exact tag or the short hash of the tree `SPEC` points to (plus `-dirty` if it has uncommitted changes outside this
folder), or the content of its `COMMIT` file.

A test that does not finish within `TIMEOUT` is recorded with status `TIMEOUT(<limit>s)`, e.g. `TIMEOUT(2400s)`, and
`make compare` shows it as `>2400`. `make all` always runs every test, also the ones that timed out before (e.g.
WKICM_4.7 with v1.0.2, IsPrincipal_6.2 with both trees). A timeout with a commit other than the baseline (`BASELINE`),
of a test that completed before on the same machine (with any commit), is also a regression: the run prints a
`REGRESSION:` warning and the case is recorded in `regressions.txt`, to be investigated.

Which runs are needed: every test has at least one run with the baseline and one with each new commit on aerox and on
babbage (a run that times out counts). So `make all` on both machines for every new commit (the baseline lines stay
in `timings.log` and are not rerun; rerun them after a Magma upgrade or on a new machine). A new test is run with
v1.0.2 and with the current commit on every machine.

Only `timings.log` and `regressions.txt` are kept under version control; the per-test logs in `logs/`
are temporary (they are in `.gitignore`) and are deleted by the next `make all`.

The CPU time (`Cputime`) of the timed part is what is logged; the wall time of each test is in `logs/`.

Do not start two runs in the same folder at the same time: they would write the same files in `logs/` (and read each
other's results). For a second concurrent run use another `LOGDIR=...` (and `LOGFILE=...` if its results should not go
into the shared log).

## Writing a test

`tests/<Name>_<g>.<n>.m`, where `<Name>` is the intrinsic or the theme, `<g>` the group number and `<n>` the number
of the test in the group. Groups: 1 PicardGroup, 2 ideal arithmetic (SmallRepresentative, ZBasisLLL, products),
3 OverOrders and GraphOverOrders, 4 WKICM, 5 ICM and IsIsomorphic, 6 IsPrincipal and UnitGroup, 7 IntermediateIdeals
and IdealsOfIndex, 8 primes, CRT, residue fields and completions, 9 Cohen-Macaulay type and units.
Take the next free number; numbers are never reused, so that the log stays meaningful.

Template:

    /*
        Name_g.n: one line saying what is timed.
        Why: what makes the example interesting (irreducible or not, signature, class group of the maximal order,
        index of the order, singular primes with their residue fields and ramification, Cohen-Macaulay type,
        numbers of overorders/classes), and the expected size of the output.
    */
    AttachSpec(spec);
    SetSeed(1);
    SetClassGroupBounds("GRH");   // if class groups are involved
    _<x>:=PolynomialRing(Integers());
    // setup that is not timed, with the fixed inputs (see below) and their checksum:
    // checksum:=function(elts) ... end function;     (copy it from an existing test)
    // data:=[ ... ];                                  (Z-bases printed with PrintSeqAlgEtQElt)
    // chk:=checksum([ A!t : t in &cat data ]);
    // printf "### INPUT %o\n",chk;
    // assert chk eq ...;
    t0:=Cputime();
        // the timed computation
    t:=Cputime(t0);
    assert ...;                   // correctness of the output
    printf "### TIMING %o %o\n",name,t;
    print "### DONE";
    quit;

`spec` and `name` are given on the command line by the Makefile. A test attaches only `spec` and uses only
intrinsics that exist in every version the suite is meant to compare (the baseline is release v1.0.2). Keep the timed
part between ~1 second and ~10 minutes on a laptop: shorter is noise, longer makes the suite slow. `make check`
verifies the required lines.

## Fixed inputs

A test must time the same objects in every tree and every Magma version. `SetSeed(1)` does not ensure this: the
random stream depends on the Magma version, and the objects built by the package (the representative `p(g)` of a
class, the generators of a class group, `ShortElement`, the order of `PrimesAbove` or of `Classes`) change with the
code that is being timed. So:

- a test does not call `Random` (`make check` fails on `Random(`);
- objects that come from a computation (ideals, primes, short elements) are stored in the test, printed with
  `PrintSeqAlgEtQElt`, and rebuilt with `Ideal(S,[ A!t : t in zb ])`; an object derived from them in the setup (a
  power, `a*I`) is taken through `ZBasis(Ideal(S,ZBasis(...)))`, the HNF basis, which is the same in every version;
- many generic elements come from the Lehmer generator `lehmer_elts` (see IdealOps_2.3);
- `A!t` reads the coordinates in the order of the components of A: when f is reducible, the test asserts
  `[ DefiningPolynomial(K) : K in Components(A) ]`;
- the test prints `### INPUT <checksum>` of its inputs (the function `checksum`) and asserts the value, so that a
  different input makes the test fail instead of timing something else;
- the randomness inside the timed code (e.g. in CoprimeRepresentative) stays, with `SetSeed(1)` at the top; when it
  dominates the running time, the test times a few seeds on fresh copies (IsPrincipal_6.2).

The header of each test says where its inputs come from. A test whose inputs change gets a new number.
