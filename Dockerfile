FROM avianalter/samap:latest

# Copy your script into the container's bin directory
COPY scripts/run_samap.py /usr/local/bin/run_samap.py
# Make sure it’s executable (optional)
RUN chmod +x /usr/local/bin/run_samap.py

WORKDIR /workspace
