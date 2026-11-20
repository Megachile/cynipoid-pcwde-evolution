# Paper build launcher

`launch.py` is the successor to task-local manuscript build wrappers. Supply a
checked JSON job configuration and run it on the login node under the project
group, for example `sg PROJECT_GROUP -c 'python scripts/build/launch.py
config.json production'`. The launcher checks the pinned repository and output
before submitting `job.sbatch` through the provenance wrapper. The Slurm job
rechecks them, stages an immutable deposit to node-local storage, and runs
`make verify`, `make manuscript`, and `make stage-check` with the repository's
unchanged targets. `fixture` mode accepts a synthetic deposit override and
tests the figure-intent gate and copy modes. Production rejects overrides.

The configuration is task-owned and must pin `repo`, `base_sha`, `sha`,
`main_git_dir`, `deposit`, `destination`, `receipt`, `prior_data`, `workspace`,
`python`, `envs`, `provenance_sbatch`, and `project_gid`. It may set `make_target`
to the documented manuscript alias. The destination must be a new child of the
checkout parent. The config and submission script are frozen and hashed by the
caller before submission. The caller retains output and receipt on project
storage; node-local staging is removed only after a successful build.
Run from a clean isolated checkout with `PYTHONDONTWRITEBYTECODE=1`;
preflight checks ignored files as well as tracked and untracked files.

For storage-bounded builds, optional `bounded_staging` config explicitly declares
`source`, `inventory` and `inventory_sha256`. The inventory is a sealed JSON list
of relative name, byte size and SHA256 for every source file. The source must be
a task-owned `staged_input_*` directory immediately next to the checkout. This
mode performs the unchanged complete `make stage-check` before rendering, checks
that the source stayed unchanged, and releases only the successful comparison
stage. It then performs the same node-local copy and two-sided hash checks as
the default route. After the source, node copy and sealed inventory all agree,
production consumes the declared task-owned source copy before rendering.
Fixture mode preserves that source for the subsequent production run. Cleanup
receipts include paths, bytes, file counts and validation gates. A failed compare,
changed source/copy/seal or unsafe path fails before deletion. Default callers
retain their sources and the original staging order. No shared launcher install
is needed. `tests/test_lifecycle.py` covers mutation and confinement failures.
