# Use black to format the codebase
python -m black `
    --line-length 100 `
    --include summer_py/**/*.py tests/**/*.py `
    .
