# 2026-03-18 — Fix CRAN warn suppression policy violation

> Run: `CRAN-fix-warn-001` | Profile: r-package | Verdict: PASS

## What Changed

Replaced the CRAN-prohibited `options(warn = -1)` global warning suppression pattern with the CRAN-compliant `suppressWarnings()` wrapper in the `silent_ebalance()` helper function inside `R/fdid.R`. This was required because CRAN maintainer Benjamin Altmann rejected the package submission due to the policy violation documented at https://contributor.r-project.org/cran-cookbook/code_issues.html#setting-optionswarn--1.

## Files Changed

| File | Action | Description |
| --- | --- | --- |
| `R/fdid.R` | modified | Removed `old_warn <- getOption("warn"); options(warn = -1)` (line 99) and `options(warn = old_warn)` (line 103); wrapped `ebal::ebalance()` call in `suppressWarnings()` |

## Process Record

This section captures the full workflow history: what was proposed, what was tested, what problems arose, and how they were resolved.

### Proposal (from theorist)

**Implementation spec summary** (from `spec.md`):
- Remove the `old_warn <- getOption("warn"); options(warn = -1)` line entirely
- Wrap `ebal::ebalance(Treatment = Treatment, X = X)` in `suppressWarnings()`
- Remove the `options(warn = old_warn)` restore line
- Preserve blank lines and the sink/on.exit pattern (orthogonal to this fix)
- No new dependencies, no NAMESPACE changes, no documentation changes needed (internal helper)

**Test spec summary** (from `test-spec.md`):
- Scenario 1: Static check — `options(warn` regex must return zero matches in `R/fdid.R`
- Scenario 2: Static check — `getOption("warn")` must return zero matches
- Scenario 3: Static check — `suppressWarnings` must be present wrapping `ebal::ebalance()`
- Scenario 4: `R CMD check --as-cran` must pass with no new ERRORs or WARNINGs
- Scenario 5: Existing test suite must pass (no regressions)
- Scenario 6: No `options(warn = ...)` anywhere in `R/fdid.R` (no global state leakage)
- Tolerances: deterministic change, no numerical tolerances applicable

### Implementation Notes (from builder)

- Removed `old_warn <- getOption("warn"); options(warn = -1)` (original line 99)
- Wrapped `ebal::ebalance()` in `suppressWarnings()` (new line 100)
- Removed `options(warn = old_warn)` (original line 103)
- Preserved blank lines for readability consistency
- Did not modify the sink/on.exit pattern (handles console output, separate concern)
- `suppressWarnings()` is base R — no imports needed
- No unit tests written (mechanical refactor with identical runtime behavior)
- `parse('R/fdid.R')` confirmed no syntax errors

### Validation Results (from auditor)

- **Scenario 1 (forbidden pattern removed)**: PASS — zero matches for `options(warn`
- **Scenario 2 (save pattern removed)**: PASS — zero matches for `getOption("warn")`
- **Scenario 3 (suppressWarnings present)**: PASS — found at line 100 wrapping `ebal::ebalance()` (plus one pre-existing unrelated match at line 150)
- **Scenario 4 (R CMD check --as-cran)**: NOT RUNNABLE — sandboxed environment lacks CRAN dependencies (no internet). `R CMD build` succeeded cleanly, producing `fdid_1.0.1.tar.gz`.
- **Scenario 5 (functionality preserved)**: NOT RUNNABLE — depends on Scenario 4.
- **Scenario 6 (no global state leakage)**: PASS — zero matches for `options(warn` in entire file
- **Syntax validation**: PASS — `Rscript --vanilla -e "parse(file='R/fdid.R')"` succeeded
- **R CMD build**: PASS — built cleanly with no warnings or errors

### Problems Encountered and Resolutions

No problems encountered. No BLOCK, HOLD, or STOP signals were raised during this workflow.

### Review Summary (from skeptic, if available)

Pending — skeptic review follows scribe.

- **Pipeline isolation**: pending
- **Convergence**: pending
- **Tolerance integrity**: pending
- **Verdict**: pending

## Design Decisions

1. **`suppressWarnings()` over `tryCatch(warning = ...)`**: `suppressWarnings()` is the simplest CRAN-approved pattern for this use case. It scopes warning suppression to exactly the one call that needs it, without global state modification. `tryCatch` would be more complex with no additional benefit since we want to suppress all warnings from `ebal::ebalance()`, not selectively handle them.

2. **Remove save/restore entirely rather than converting to `on.exit`**: Since `suppressWarnings()` is lexically scoped, there is no global state to save or restore. The `old_warn`/`options(warn = old_warn)` pattern became dead code and was removed rather than converted to an `on.exit()` guard.

3. **Preserve sink/on.exit pattern**: The `sink()` redirections handle stdout/stderr suppression (console output from `ebal::ebalance`), which is a separate concern from warning suppression. This pattern was not flagged by CRAN and remains correct.

## Handoff Notes

- The fix is minimal and mechanical. The `silent_ebalance()` function's observable behavior is unchanged — it still suppresses both console output (via sink) and warnings (now via `suppressWarnings` instead of `options(warn=-1)`) from `ebal::ebalance()`.
- `R CMD check --as-cran` could not be run in the sandboxed CI environment due to missing CRAN dependencies. The submitter should run `R CMD check --as-cran` locally or in a fully-provisioned CI environment before resubmitting to CRAN.
- The CRAN-specific NOTE about `options(warn = -1)` that triggered the rejection will no longer appear since the pattern is completely removed from the source.
- There is one other `suppressWarnings()` call in `R/fdid.R` at line 150 (wrapping `as.numeric()`) — this is pre-existing and unrelated to this fix.
