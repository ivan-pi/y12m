# Fortran AGENT Guidelines

## Licensing

- All contributions must be compatible with **GPL-2.0-only**.
- Every new source file (Fortran, C, CMake, shell) must start with an SPDX
  license header on the first line:
  - Fortran / shell: `! SPDX-License-Identifier: GPL-2.0-only`
  - C: `// SPDX-License-Identifier: GPL-2.0-only`
  - CMake / script: `# SPDX-License-Identifier: GPL-2.0-only`

## Attribution

When an AI tool contributes code, include an attribution comment near the top
of each new file (after the SPDX header):

```
! Assisted-by: AGENT_NAME:MODEL_VERSION
```

For example: `! Assisted-by: GitHub Copilot:claude-sonnet-4.5`

Do **not** add any "Signed-off-by" tags; only the human reviewer can certify
authorship.

## Fortran Style preferences

- Separate statements on different lines; do not use semi-colons (;).
- Do not use implicit casting; use the `real()` and `int()` intrinsic functions; do not use the legacy forms like `aint` or `dble`.
- Use `real([kind=..])` to declare real variables. DON'T use the legacy `double precision` form, unless you are editing an older file that already uses it.
- For real variables use a `sp` (single precision), `dp` (double precision) or `wp` (working precision) kind specifier; an exception is in files using the C interop features.
- Do not bother aligning the type declarations and variables with the `::` token; just use one space on the left and right.
- Use 3 spaces for indentation by default or adapt to the existing tab size used in a given file.
- Always add `implicit none` in interface blocks for external or `bind(C)` procedures; you can skip this rule when a procedure only has a few parameters and the chance of missing one is low.
- Always use the `only` clause on module `use` statements.
- Always add `intrinsic` in case of the `iso_c_binding` and `iso_fortran_env` modules.
- When generating new procedures, always make use of the `intent(..)` attributes.
- For procedures that use pass by `value`, you can omit the `intent`.
- Modules should use `private` attribute and use explicit `public` statements to export only the minimum necessary.
- Any new code should use free-form Fortran and the `.f90`/`.F90` file extensions.

## Repository conventions

- **Preserve the legacy/fixed-form boundary.** Do not rewrite files under
  `src/legacy/*.f` to free-form unless explicitly requested; keep new code in
  separate `.f90`/`.F90` files.
- **Edit templates, not generated sources.** The files `y12ma.f`, `y12mb.f`,
  `y12mc.f`, `y12md.f`, `y12mg.f`, and `y12mh.f` in `src/legacy/` are
  generated from the Fypp templates in `templates/`.  Never edit these `.f`
  files directly: change the corresponding `.fpp` template, regenerate with
  `make -C templates` (requires `pip install fypp`), and commit the template
  together with the regenerated source.  CI fails if the two diverge.  Only
  the iterative-refinement pair is maintained by hand: `y12mfe.f` (single
  precision, fixed-form) and `src/y12mff.F90` (double precision, free-form;
  excluded from the template pattern because its residual accumulator —
  quad precision or double-double — cannot be expressed as a simple kind
  substitution).
- **Prefer module wrappers over changing legacy entry points.** Put modern
  interface improvements in `src/y12m.F90` or helper modules; only edit a
  legacy routine when the bug genuinely lives there.
- **Keep single- and double-precision support aligned.** If a new helper,
  test, or interface is meant to be generic, add both `sp` and `dp` variants
  unless there is a documented reason not to.
- **Reuse shared helpers.** Extend `test/y12m_test_helpers.f90` and
  `examples/y12m_example_helpers.f90` before adding ad-hoc local copies of
  the same functionality.
- **Use absolute paths in CMake tests** that refer to source-tree data files
  (e.g. `.mtx` inputs), so tests work regardless of the build directory.

## Documentation

- When adding or changing public APIs, update `README.md`, `docs/API.md`, and
  `examples/README.md` (API usage table) in the same change.
- Keep comments high-value and technical: prefer notes that explain solver
  usage constraints, API ordering requirements, or numerical assumptions rather
  than line-by-line restatements of obvious code.

## Testing instructions

- Add or update tests for the code you change, even if nobody asked.
- Any bug fix must include a regression test that fails before the fix and
  passes after.
- Build and run the suite with:
  ```sh
  cmake -B build && cmake --build build && ctest --test-dir build
  ```
- Prefer self-checking example programs (manufactured solutions, known exact
  solutions, or explicit error metrics) over programs that merely run without
  crashing.
