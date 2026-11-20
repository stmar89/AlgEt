# Timing tests

A suite of timing tests for the package AlgEtQ (the files listed in `spec`), one test per file, run in parallel on a
limited number of cores, with the timings logged per machine, Magma version and commit. It replaces the file of slow
tests, now kept in `dev/old/slow_tests_AlgEtQ.m`. AlgEtNF and the packages of `specMod` (AlgEtQMod) and `specMtrx`
(AlgEtQMtrx) are not covered. The baseline is commit d65a92f (until 2026-11-20 it was release v1.0.2): see History
below.

## Running

From this folder:

    make all                                   # all the tests, 4 at a time, against ../../spec
    make all JOBS=2 MACHINE="laptop aerox"     # fewer parallel jobs, a name for the machine in the log
    make all SPEC=/path/to/other/tree/spec     # the same tests against another checkout (e.g. a `git archive`
                                               # of a release: put the commit hash in a file COMMIT next to its spec)
    make run TESTS="PicardGroup_1.7 WKICM_4.*" # a subset
    make summary                               # the table of the last run again
    make compare A=d65a92f B=feature MACHINE="laplace OvGU"  # medians from timings.log, side by side
    make list | make check | make clean

Variables: `JOBS` (parallel tests, default 4: keep it well below the number of cores, so that the timings do not
interfere), `SPEC` (default `../../spec`), `MACHINE` (free text, default `hostname`: eg. "laptop aerox",
"babbage" and "laplace OvGU"), `LABEL` (a free comment stored in the log), `TIMEOUT` (seconds per test, default 2400),
`TESTS`, `BASELINE` (the reference commit, default d65a92f, see below).

Each run appends one line per test to `timings.log` (tab separated: date, machine, Magma version, commit, label,
test, CPU seconds, status). The date is in CET (UTC+1 all year, no summer time) on every machine. The commit is the
exact tag or the short hash of the tree `SPEC` points to (plus `-dirty` if it has uncommitted changes outside this
folder), or the content of its `COMMIT` file.

A test that does not finish within `TIMEOUT` is recorded with status `TIMEOUT(<limit>s)`, e.g. `TIMEOUT(2400s)`, and
`make compare` shows it as `>2400`. `make all` always runs every test, also the ones that timed out before (e.g.
WKICM_4.7 with v1.0.2). A timeout with a commit other than the baseline (`BASELINE`),
of a test that completed before on the same machine (with any commit), is also a regression: the run prints a
`REGRESSION:` warning and the case is recorded in `regressions.txt`, to be investigated.

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
Take the next free number; numbers are never reused, so that the log stays meaningful. A retired test is moved to
`tests/old/` (only the files directly in `tests/` are run).

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
intrinsics that exist in every version the suite is meant to compare (the baseline is commit d65a92f). Keep the timed
part between 2 and 10 minutes with the baseline, on a laptop or on laplace (they run at about the same speed):
shorter is noise, longer makes the suite slow. To get there, a test can time a family of examples (several algebras
or orders, with different arithmetic: degrees, signatures, splitting of the singular primes, class groups) rather than
one example; repeating the same computation on fresh copies of the algebra is the last resort. `make check` verifies
the required lines.

## Fixed inputs

A test must time the same objects in every tree and every Magma version. `SetSeed(1)` does not ensure this: the
random stream depends on the Magma version, and the objects built by the package (the representative `p(g)` of a
class, the generators of a class group, `ShortElement`, the order of `PrimesAbove` or of `Classes`) change with the
code that is being timed. So:

- a test does not call `Random` (`make check` fails on `Random(`);
- objects that come from a computation (ideals, primes, short elements) are stored in the test, printed with
  `PrintSeqAlgEtQElt`, and rebuilt with `Ideal(S,[ A!t : t in zb ])`; an object derived from them in the setup (a
  power, `a*I`) is taken through `ZBasis(Ideal(S,ZBasis(...)))`, the HNF basis, which is the same in every version;
- many generic elements come from the Lehmer generator `lehmer_elts` (see IdealOps_2.6), and generic invertible
  ideals are the first N invertible ideals generated by two of its elements (`lehmer_inv_ideals`, see
  PicardGroup_1.10); ideals with a known shape are built from integers and the primitive element F, e.g. the primes
  of degree 1 lE + (F-c)E for l not dividing the index (IsPrincipal_6.3), and their powers and products are taken
  through their HNF bases;
- objects that the code returns as a set in an order that may change (the primes above the conductor, the overorders
  of a given Cohen-Macaulay type) are sorted by a key that does not depend on the code, e.g. the minimal integer, the
  index and the HNF basis (CRT_8.3), or used only as a set (WKICM_4.16);
- `A!t` reads the coordinates in the order of the components of A: when f is reducible, the test asserts
  `[ DefiningPolynomial(K) : K in Components(A) ]`;
- the test prints `### INPUT <checksum>` of its inputs (the function `checksum`) and asserts the value, so that a
  different input makes the test fail instead of timing something else;
- the randomness inside the timed code (e.g. in CoprimeRepresentative) stays, with `SetSeed(1)` at the top; when it
  dominates the running time, the test times a few seeds on fresh copies (PicardGroup_1.7);
- `Order` returns the order already known to the algebra (the attribute `KnownOrders`), with all its attributes, so
  a fresh copy of an order needs a fresh copy of the algebra (`EtaleAlgebra(f)` again). The graph tests use this to
  time GraphOverOrders when the overorders are known but not their minimal inclusions, as after LoadWKICM: the
  overorders are computed in one copy and rebuilt from their Z-bases in another (`copy_with_overorders`, see
  GraphOverOrders_3.14).

The header of each test says where its inputs come from. A test whose inputs change gets a new number.

## History

On 2026-11-20 the baseline moved from release v1.0.2 to commit d65a92f, and the 35 tests that took less than 20
seconds with 9362a5e on laplace were retired: each was replaced by a test with a new number, in the same spirit but
with larger examples or families of examples, so that it takes 2 to 10 minutes. The retired tests are kept in
`tests/old/`, which `make` does not run, and their lines stay in `timings.log` (and in `regressions.txt`). The 11
tests that took longer were kept: OverOrders_3.7, WKICM_4.7, WKICM_4.8, WKICM_4.9, WKICM_4.10, WKICM_4.11, WKICM_4.14,
WkAbstract_5.4, ICMAbstract_5.5, IntermediateIdeals_7.1 and IdealsOfIndex_7.2.

| retired | new | | retired | new | | retired | new |
|---|---|---|---|---|---|---|---|
| PicardGroup_1.1 | PicardGroup_1.7 | | OverOrders_3.8 | OverOrders_3.13 | | WKICM_4.12 | WKICM_4.21 |
| PicardGroup_1.2 | PicardGroup_1.8 | | GraphOverOrders_3.4 | GraphOverOrders_3.14 | | WKICM_4.13 | WKICM_4.22 |
| PicardGroup_1.3 | PicardGroup_1.9 | | GraphOverOrders_3.5 | GraphOverOrders_3.15 | | ICM_5.1 | ICM_5.6 |
| PicardGroup_1.4 | PicardGroup_1.10 | | GraphOverOrders_3.6 | GraphOverOrders_3.16 | | ICM_5.2 | ICM_5.7 |
| PicardGroup_1.5 | PicardGroup_1.11 | | GraphOverOrders_3.9 | GraphOverOrders_3.17 | | ICM_5.3 | ICM_5.8 |
| PicardGroup_1.6 | PicardGroup_1.12 | | WKICM_4.1 | WKICM_4.15 | | IsPrincipal_6.1 | IsPrincipal_6.3 |
| SmallRepresentative_2.1 | SmallRepresentative_2.4 | | WKICM_4.2 | WKICM_4.16 | | IsPrincipal_6.2 | IsPrincipal_6.4 |
| ZBasisLLL_2.2 | ZBasisLLL_2.5 | | WKICM_4.3 | WKICM_4.17 | | CRT_8.1 | CRT_8.3 |
| IdealOps_2.3 | IdealOps_2.6 | | WKICM_4.4 | WKICM_4.18 | | Completion_8.2 | Completion_8.4 |
| OverOrders_3.1 | OverOrders_3.10 | | WKICM_4.5 | WKICM_4.19 | | CohenMacaulay_9.1 | CohenMacaulay_9.4 |
| OverOrders_3.2 | OverOrders_3.11 | | WKICM_4.6 | WKICM_4.20 | | CohenMacaulay_9.2 | CohenMacaulay_9.5 |
| OverOrders_3.3 | OverOrders_3.12 | | | | | TotRealPos_9.3 | TotRealPos_9.6 |

The header of each new test says which test it replaces.
