<!--
BASE BRANCH: pull requests go to `development`, NOT `main`.
`main` only receives merges from `development` at release time, so a PR based on
`main` cannot be merged as-is. GitHub preselects `main` because it is this
repository's default branch -- if the base shown above is `main`, click "Edit"
beside the title and change it to `development`.

See CONTRIBUTING.md, including how to retarget and rebase a PR already opened
against `main`.
-->

## What this changes


## Checklist

- [ ] the base branch is `development`;
- [ ] new or changed behaviour is covered by a `testthat` test;
- [ ] `devtools::document()` re-run, using the roxygen2 version pinned in `DESCRIPTION`;
- [ ] `NEWS.md` bullet added under the *development version* headings (not a released version's);
- [ ] the development version suffix in `DESCRIPTION` is bumped;
- [ ] `R CMD check --as-cran` passes.
