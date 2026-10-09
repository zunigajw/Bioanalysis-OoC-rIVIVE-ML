FROM python:3.12-slim

# make runs the pipeline; libxrender1/libxext6 are needed by the RDKit wheel.
RUN apt-get update && \
    apt-get install -y --no-install-recommends make ca-certificates libxrender1 libxext6 && \
    rm -rf /var/lib/apt/lists/*

# Python packages, pinned for reproducibility.
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Show each script's progress as it runs.
ENV PYTHONUNBUFFERED=1
WORKDIR /work
