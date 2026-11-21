# Python-only image of the analytical package (NO MOOSE inside; MOOSE stays an external, optional solver).
FROM python:3.12.9-slim
ENV PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
WORKDIR /opt/yield_surface_discovery
COPY requirements-lock.txt .
RUN pip install --no-cache-dir -r requirements-lock.txt
COPY . .
RUN python scripts/validate_install.py && python -m pytest -q
CMD ["python", "scripts/run_minimal_example.py", "--check", "--output", "/tmp/minimal_example"]
