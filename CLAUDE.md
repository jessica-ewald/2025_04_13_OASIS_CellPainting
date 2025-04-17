# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build/Run Commands
- Install dependencies: `uv sync`
- Run Snakemake workflow: `cd 01_snakemake && snakemake --cores 4`
- Run single Snakemake rule: `cd 01_snakemake && snakemake --cores 1 rule_name`
- Run specific output file: `cd 01_snakemake && snakemake --cores 1 path/to/target_file`
- Create workflow DAG: `cd 01_snakemake && snakemake --dag | dot -Tpdf > workflow.pdf`
- Lint Python code: `ruff check .`
- Format Python code: `ruff format .`

## Code Style Guidelines
- Use Python 3.13+ features and syntax
- Formatting: Follow ruff formatter conventions (similar to black)
- Imports: Group imports (stdlib, third-party, local) and sort alphabetically
- Naming: snake_case for functions/variables, PascalCase for classes
- Type hints: Use proper type annotations for function parameters and returns
- Error handling: Use specific exceptions with descriptive messages
- Documentation: Add docstrings for modules, classes, and non-trivial functions
- When modifying R files, maintain existing formatting style
- Use pandas/polars for data manipulation and numpy for numerical operations
- Follow Snakemake best practices for rule definitions