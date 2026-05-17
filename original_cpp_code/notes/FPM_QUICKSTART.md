# FPM Quick Reference for Jacchia-Roberts

## Building

```bash
# Build everything (library + executables)
pixi run fpm build

# Build with optimization (release profile)
pixi run fpm build --profile release

# Build with debug symbols
pixi run fpm build --profile debug
```

## Running

```bash
# Run the test program
pixi run fpm run test_jacchia_roberts

# Run with specific compiler
pixi run fpm run --compiler gfortran

# Run all executables
pixi run fpm run
```

## Testing

```bash
# Run tests (when you add test/ directory)
pixi run fpm test

# Run specific test
pixi run fpm test <test_name>
```

## Development

```bash
# Enter a shell with environment set up
pixi shell
fpm build
fpm run

# Clean build artifacts
rm -rf build/
```

## Custom Compiler Flags

Edit `fpm.toml` to add custom flags:

```toml
[build]
auto-executables = true
external-modules = []

[fortran]
implicit-typing = false
source-form = "free"

[profile.release]
compiler = "gfortran"
flags = "-O3 -march=native"

[profile.debug]
compiler = "gfortran"
flags = "-g -Wall -Wextra -fcheck=all"
```

## Package Information

```bash
# Show package info
fpm --version
fpm list

# Show dependencies
fpm list --list
```

## Integration

Use as a dependency in your `fpm.toml`:

```toml
[dependencies]
jacchia-roberts = { path = "/path/to/jacchia-roberts-fortran" }
```

Or from a git repository:

```toml
[dependencies]
jacchia-roberts = { git = "https://github.com/user/jacchia-roberts-fortran.git" }
```

## Files Structure

```
jacchia-roberts-fortran/
├── fpm.toml                    # Package manifest
├── src/
│   └── jacchia_roberts_module.f90  # Library source
└── test_jacchia_roberts.f90    # Test executable
```

## Notes

- FPM automatically handles module dependencies
- Build artifacts go to `build/` directory (gitignored)
- Executables are placed in `build/gfortran_*/app/`
- Library is compiled to `build/gfortran_*/library/`
