#!/bin/bash
set -e
# Install pyphylon from bind-mounted source (fast, pure Python)
cd /pipeline && pip install -e . -q 2>/dev/null
# Run the requested command
exec "$@"
