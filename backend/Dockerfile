FROM python:3.11-slim

WORKDIR /app

# System deps required by RDKit
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements-prod.txt .
RUN pip install --no-cache-dir -r requirements-prod.txt

COPY . .

# Vina binary must be executable
RUN chmod +x bin/vina_linux_x86_64

# SQLite data directory — writable by HF Spaces non-root user
RUN mkdir -p data && chmod 777 data

# HF Spaces requires port 7860
EXPOSE 7860

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "7860"]
