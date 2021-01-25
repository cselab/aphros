# Contribution guide

## Code formatting

### C++

Run clang-format for C/C++ files found recursively in current directory

```
ap.format
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

