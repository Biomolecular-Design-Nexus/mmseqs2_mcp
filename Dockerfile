FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-venv \
    wget tar \
    && rm -rf /var/lib/apt/lists/*

# Symlink python3 -> python so scripts work uniformly
RUN ln -sf /usr/bin/python3 /usr/bin/python

WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --no-cache-dir -U cryptography certifi

# Download and install MMseqs2 GPU binary
RUN wget -q https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz \
    && tar xzf mmseqs-linux-gpu.tar.gz \
    && cp mmseqs/bin/mmseqs /usr/local/bin/ \
    && rm -rf mmseqs mmseqs-linux-gpu.tar.gz \
    && mmseqs version

# Copy application source
COPY src/ ./src/
COPY examples/ ./examples/

# Create working directories
RUN mkdir -p tmp/inputs tmp/outputs

ENV PYTHONPATH=/app

# Default database path (mount your database volume here)
ENV MMSEQS2_DB_PATH=/data/uniref100/uniref100.fasta.db_padded

CMD ["python", "src/server.py"]
