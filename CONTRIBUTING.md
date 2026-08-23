# Contributing to LandR

We welcome contributions to `LandR`, both as bug reports and as package enhancements.

## Report a bug

1. Use the [issue tracker](https://github.com/PredictiveEcology/LandR/issues) to report a bug.

2. Please include a minimal [reproducible example](https://stackoverflow.com/q/5963269/1380598) that triggers the bug.

3. Please include the output of `devtools::session_info()`.

## Branching model

This repository uses the [Git Flow](https://nvie.com/posts/a-successful-git-branching-model/) branching model (the [`git flow`](https://github.com/petervanderdoes/gitflow-avh) extension is useful for this).
Changes flow in one direction only:

```
feature / bugfix / test branch  -->  development  -->  main
```

- [`main`](https://github.com/PredictiveEcology/LandR/tree/main) holds the code of the latest release.
  It only ever receives merges from `development` at release time, and is **not** a destination for contributions.
- [`development`](https://github.com/PredictiveEcology/LandR/tree/development) holds the latest contributions and everything else queued for the next release.
- Feature, bugfix, and test branches are cut from `development` and merged back into `development`.

`main` is this repository's default branch, so GitHub preselects it as the base when you open a pull request.
**Change the base branch to `development` before you submit.**

Other long-lived branches you may see (e.g. `dev-stable`, `LandWeb`) are project-specific pins; don't target them unless a maintainer asks you to.

## Submit an enhancement

Branch from an up-to-date `development`:

```bash
## once, if you are working from a fork:
git remote add upstream https://github.com/PredictiveEcology/LandR.git

git fetch upstream
git switch --create my-feature upstream/development
```

Then send a [pull request](https://docs.github.com/articles/using-pull-requests/) with `development` as the destination branch.

### Before you open the pull request

- new or changed behaviour is covered by a `testthat` test in `tests/testthat/test-<file>.R`, matching the `R/<file>.R` it exercises;
- roxygen documentation is updated and regenerated with `devtools::document()`.
  The roxygen2 version is pinned in `DESCRIPTION` (`Config/roxygen2/version`) -- please use that version, so unrelated files in `man/` don't churn;
- `NEWS.md` has a bullet under the *development version* headings at the top of the file (`## Bug fixes`, `## New features`, ...), **not** under a released version's heading;
- the development version suffix in `DESCRIPTION` is bumped (e.g. `1.2.0.9007` -> `1.2.0.9008`);
- `R CMD check --as-cran` passes (`devtools::check(cran = TRUE)`); it also runs automatically on the pull request;
- the diff is limited to your change -- please don't reformat surrounding or otherwise untouched code, since that buries the substance of the change.

### I already opened a pull request against `main`

There's no need to close it and start over:

1. On the pull request page, click **Edit** beside the title and change the base branch to `development`.

2. Rebase your branch onto `development` so it doesn't carry commits that exist only on `main`:

   ```bash
   git fetch upstream
   git rebase --onto upstream/development upstream/main my-feature
   git push --force-with-lease
   ```

   Note the `--onto` form: a plain `git rebase upstream/development` would try to replay every commit that `main` has and `development` doesn't, which conflicts heavily.

3. Re-check the items above against `development` -- in particular, `NEWS.md` bullets often need to move, because the released-version sections on `main` don't line up with the development-version sections on `development`.

We'll try to review your pull request and provide feedback / merge improvements as quickly as possible.

Thank you!
