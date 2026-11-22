FROM public.ecr.aws/docker/library/gcc:13-bookworm@sha256:930f2ebe239275fa67226654cb79273ea34eee672ae61c8a39f689c37fb7ac5c

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pytest tmux asciinema ca-certificates git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY task_file/input_data/ /app/input_data/
COPY task_file/scripts/baseline.cpp /app/scripts/baseline.cpp
RUN cp /app/scripts/baseline.cpp /app/solve.cpp
