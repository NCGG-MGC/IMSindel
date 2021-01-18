FROM ruby:3.0.0-slim

RUN set -eux; apt-get update && \
    apt-get install -y  --no-install-recommends samtools curl mafft && \
    gem install bio && \
    cd /opt && \
    curl -OL http://faculty.virginia.edu/wrpearson/fasta/CURRENT/fasa36-linux64.tar.gz && \
    tar xf fasa36-linux64.tar.gz && \
    mkdir /opt/IMSindel && \
    rm -rf /var/lib/apt/lists/*

ENV PATH $PATH:/opt/fasta-36.3.8h/bin
COPY . /opt/IMSindel/
ENTRYPOINT ["/opt/IMSindel/bin/imsindel", "--temp", "/dev/shm", "--thread", "1"]
CMD ["--help"]
