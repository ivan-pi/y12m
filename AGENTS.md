# Fortran AGENT Guidelines

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

## Testing instructions
- Add or update tests for the code you change, even if nobody asked.
