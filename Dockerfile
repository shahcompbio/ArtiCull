# syntax=docker/dockerfile:1
#
# Lean runtime image for ArtiCull.
#
#   docker build -t articull:latest .
#   docker run --rm --shm-size=2g \
#       -v "$PWD:/data" \
#       -v /path/to/resources:/opt/articull/resources:ro \
#       articull:latest classify \
#           /data/sample.maf /data/out/sample /opt/articull/models/preprint_model /data/sample.bam -j 8
#
# Python is capped at 3.11: scikit-learn 1.2.2 -- required to unpickle the
# shipped model -- publishes no wheels for 3.12+.

# ---------------------------------------------------------------- builder ---
FROM python:3.11-slim-bookworm AS builder

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Deliberately smaller than requirements.yml. Dropped:
#   R, R-essentials        - only used by train_preprocessing, via
#                            scripts/misc/extract_cell2clone.R
#   matplotlib, seaborn    - only imported by articull.train.generate_labels
#                            and articull.train.preprocessing
#   ucsc-bigwigtobedgraph  - builds the mappability resources; not a Python dep
#   bedtools               - a system package, installed in the runtime stage
#
# scikit-learn is pinned exactly: models/preprint_model/model.pkl was pickled
# with 1.2.2, and unpickling under a different minor version is unsupported.
# numpy is capped below 2.0 because the scikit-learn 1.2.2 wheels are built
# against the numpy 1.x C ABI.
#
# Everything lands in a self-contained venv so the runtime stage can take it as
# a single COPY and inherit none of the build-time cruft.
RUN python -m venv /opt/venv \
 && /opt/venv/bin/pip install --no-cache-dir \
        'numpy>=1.24,<2' \
        'scipy>=1.10,<1.12' \
        'pandas>=1.5,<2.0' \
        'scikit-learn==1.2.2' \
        'pysam>=0.21' \
        'pandarallel>=1.6' \
        'joblib>=1.2' \
        'psutil>=5.9' \
 && /opt/venv/bin/pip uninstall -y pip setuptools wheel \
    # ~90 MB of test suites that nothing in articull imports.
 && find /opt/venv -type d -name tests -prune -exec rm -rf {} + \
 && find /opt/venv -type d -name __pycache__ -path '*/tests/*' -prune -exec rm -rf {} +

# ---------------------------------------------------------------- runtime ---
FROM python:3.11-slim-bookworm

ARG TARGETARCH

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH=/opt/venv/bin:$PATH

# bedtools is the only external binary the classify path shells out to
# (`bedtools sort` / `bedtools map` in articull/classify/extract_features.py).
# curl + ca-certificates are for articull/setup_mappability_track.bash.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        bedtools \
        curl \
        ca-certificates \
 && rm -rf /var/lib/apt/lists/* \
    # The venv is the only interpreter environment used at runtime.
 && python3 -m pip uninstall -y pip setuptools wheel 2>/dev/null || true

# bigWigToBedGraph is needed only to *build* the mappability resources, never to
# classify. UCSC ships linux binaries for x86_64 only; on arm64, build the
# bedGraphs elsewhere and mount them.
RUN if [ "$TARGETARCH" = "amd64" ]; then \
        curl -fsSL -o /usr/local/bin/bigWigToBedGraph \
            https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph \
     && chmod +x /usr/local/bin/bigWigToBedGraph; \
    else \
        echo "Skipping bigWigToBedGraph: no UCSC build for linux/$TARGETARCH"; \
    fi

COPY --from=builder /opt/venv /opt/venv

# The package is not pip-installed: articull/__init__.py has neither a docstring
# nor __version__, so the flit backend in pyproject.toml cannot build it. Adding
# both would let this become `pip install .`.
WORKDIR /opt/articull
COPY articull/ ./articull/
COPY models/preprint_model/ ./models/preprint_model/
ENV PYTHONPATH=/opt/articull

# articull resolves its default --resources_dir to <parent of package>/resources,
# i.e. /opt/articull/resources, so mounting the mappability tracks there means
# callers never have to pass --resources_dir.
VOLUME ["/opt/articull/resources"]

WORKDIR /data
ENTRYPOINT ["python", "-m", "articull"]
CMD ["--help"]
