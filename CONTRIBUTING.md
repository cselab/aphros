# Contribution guide

## Code formatting

### C++

Format tracked C/C++ files with the pinned version of clang-format
(requires [uv](https://docs.astral.sh/uv/))

```
tools/format
```

Check without changing files (as in CI)

```
tools/format --check
```

### Python

Format all python files found recursively in current directory

```
yapf -ir .
```

## Add copyright notice

Add copyright notice to C/C++ source files found recursively in current
directory (if `copyright` is not found in first 10 lines of the file)

```
ap.applycopyright $(ap.findsource)
```

