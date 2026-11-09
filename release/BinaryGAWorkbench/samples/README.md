# Sample selection data

The `.aop` files in this directory are small square setup-cost matrices used by
both applications and by `--demo`.

## Format

```
n
c00 c01 ... c0(n-1)
c10 c11 ... c1(n-1)
...
c(n-1)0 ... c(n-1)(n-1)
```

* The first non-comment line is the integer matrix size `n` (`2 <= n <= 200`).
* The next `n` lines each contain `n` numbers. Any run of spaces or tabs is a
  valid separator, and a comma is accepted as well.
* Blank lines and lines beginning with `#` are ignored.
* Numbers are parsed with the invariant culture (`12.5`, not `12,5`).
* `c[i][j]` is the setup cost of assigning job `i` to machine `j`. Every cell is
  meaningful, including the diagonal: if the diagonal is zero then the identity
  selection is free, so the sample diagonals deliberately use non-zero values
  to keep the optimisation interesting.

## Files

| File | Size | Purpose |
| --- | --- | --- |
| `matrix5.aop` | 5 x 5 | Default demo matrix; small enough for an exhaustive optimum. |
| `matrix8.aop` | 8 x 8 | Larger matrix for longer runs and screenshots. |

Both demo programs look for `samples\matrix5.aop` next to the executable,
then in each parent directory, and finally in the current working directory.
Use `--sample <path>` to load a different file.

