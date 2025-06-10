FROM docker.io/avianalter/samap:latest

# Copy scripts into the container's bin directory
COPY scripts/run_samap.py /usr/local/bin/run_samap.py
COPY scripts/visualize_samap.py /usr/local/bin/visualize_samap.py
# Make sure it’s executable
RUN chmod +x /usr/local/bin/run_samap.py
RUN chmod +x /usr/local/bin/visualize_samap.py

# Load the custom patch to fix the analysis module
COPY patches/analysis.py /root/miniconda/lib/python3.8/site-packages/samap/analysis.py

WORKDIR /workspace
