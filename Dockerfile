FROM ruby:2.7-slim

RUN set -eux; apt-get update && \
    apt-get install -y  --no-install-recommends samtools mafft curl && \
    gem install bio && \
    cd /opt && \
    curl -OL http://faculty.virginia.edu/wrpearson/fasta/CURRENT/fasa36-linux64.tar.gz && \
    tar xf fasa36-linux64.tar.gz && \
    mkdir /opt/IMSindel && \
    rm -rf /var/lib/apt/lists/*

COPY . /opt/IMSindel/
ENV PATH $PATH:/opt/fasta-36.3.8h/bin
ENTRYPOINT ["/opt/IMSindel/bin/imsindel", "--temp", "/dev/shm", "--thread", "2"]
CMD ["--help"]
