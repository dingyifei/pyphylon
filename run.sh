#!/usr/bin/env bash
set -euo pipefail
CMD="${1:-help}"

HTML_DIR="output/html"
BUILD_STAMP=".docker-build-stamp"

# Auto-build image if Dockerfile or requirements changed since last build
_ensure_built() {
  if [ ! -f "$BUILD_STAMP" ] \
     || [ Dockerfile -nt "$BUILD_STAMP" ] \
     || [ requirements.txt -nt "$BUILD_STAMP" ]; then
    echo "=== Building Docker image ==="
    docker compose build
    touch "$BUILD_STAMP"
  fi
}

case "$CMD" in
  edit)      _ensure_built; docker compose up marimo ;;
  app)       _ensure_built; docker compose up marimo-app ;;
  pipeline)  _ensure_built; shift; docker compose run --rm pipeline snakemake --cores "${1:-all}" --sdm apptainer ;;
  export)
    _ensure_built
    mkdir -p "$HTML_DIR"
    for nb in notebooks/*.py; do
      name=$(basename "$nb" .py)
      echo "Exporting $name ..."
      docker compose run --rm pipeline marimo export html "$nb" -o "$HTML_DIR/${name}.html"
    done
    echo "HTML reports written to $HTML_DIR/"
    ;;
  run)
    _ensure_built
    CORES="${2:-all}"
    echo "=== Running snakemake pipeline (cores=$CORES) ==="
    docker compose run --rm pipeline snakemake --cores "$CORES" --sdm apptainer
    echo ""
    echo "=== Exporting notebooks to HTML ==="
    mkdir -p "$HTML_DIR"
    for nb in notebooks/*.py; do
      name=$(basename "$nb" .py)
      echo "Exporting $name ..."
      docker compose run --rm pipeline marimo export html "$nb" -o "$HTML_DIR/${name}.html"
    done
    echo "Done. HTML reports in $HTML_DIR/"
    ;;
  shell)     _ensure_built; docker compose run --rm pipeline bash ;;
  build)     docker compose build; touch "$BUILD_STAMP" ;;
  down)      docker compose down ;;
  *)
    echo "Usage: ./run.sh <command>"
    echo ""
    echo "Commands:"
    echo "  run        Full pipeline: snakemake + export HTML (pass cores as arg, default: all)"
    echo "  pipeline   Run snakemake pipeline only  (pass cores as arg, default: all)"
    echo "  export     Export all notebooks to static HTML"
    echo "  edit       Start marimo editor          (localhost:2718)"
    echo "  app        Start marimo app viewer      (localhost:2719)"
    echo "  shell      Interactive shell in container"
    echo "  build      Build/rebuild Docker image"
    echo "  down       Stop all services"
    ;;
esac
