# syntax=docker/dockerfile:1

FROM python:3.11-bookworm

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1
# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1

WORKDIR /app


# Create a non-privileged user that the app will run under.
# See https://docs.docker.com/go/dockerfile-user-best-practices/
ARG UID=10001
RUN adduser \
    --disabled-password \
    --gecos "" \
    --home "/home/peptdeepuser" \
    --shell "/sbin/nologin" \
    --uid "${UID}" \
    peptdeepuser


COPY requirements requirements
# restrict to CPU-only PyTorch wheels (to reduce image size), part 1
RUN sed -i '/^torch/i\--extra-index-url https://download.pytorch.org/whl/cpu' requirements/requirements.txt


RUN pip install --no-cache-dir  -r requirements/requirements.txt
RUN pip install --no-cache-dir  -r requirements/requirements_gui.txt

# restrict to CPU-only PyTorch wheels, part 2
RUN sed -i '/^--extra-index-url https:\/\/download.pytorch.org\/whl\/cpu$/d' requirements/requirements.txt
RUN sed -i '/^torch/d' requirements/requirements.txt

COPY peptdeep peptdeep
COPY MANIFEST.in MANIFEST.in
COPY LICENSE.txt LICENSE.txt
COPY README.md README.md
COPY pyproject.toml pyproject.toml

RUN pip install --no-cache-dir ".[stable,gui-stable]"

ENV PORT=8501
EXPOSE 8501

USER peptdeepuser

CMD ["/usr/local/bin/peptdeep", "gui", "--port", "8501"]

# build & run:
# docker build --progress=plain -t peptdeep .
# DATA_FOLDER=/path/to/local/data
# docker run -p 8501:8501 -v $DATA_FOLDER:/app/data/ -t peptdeep
