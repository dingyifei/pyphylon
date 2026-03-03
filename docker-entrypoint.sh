#!/bin/bash
set -e
# Install pyphylon from bind-mounted source (fast, pure Python)
cd /pipeline && pip install -e . -q 2>/dev/null
# Symlink project marimo.toml so marimo finds it via XDG config path
mkdir -p /root/.config/marimo
ln -sf /pipeline/marimo.toml /root/.config/marimo/marimo.toml
# Run the requested command
exec "$@"
