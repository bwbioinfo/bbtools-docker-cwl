# Use the official Rocky Linux base image
FROM rockylinux/rockylinux:latest

# Set the maintainer label
LABEL maintainer="your-email@example.com"

# Install necessary dependencies and tools
RUN yum update -y && \
    yum install -y java-1.8.0-openjdk wget && \
    yum clean all

# Set the working directory
WORKDIR /opt

# Download and install BBTools
RUN wget https://sourceforge.net/projects/bbmap/files/latest/download/BBMap_38.86.tar.gz && \
    tar -xzf BBMap_38.86.tar.gz && \
    rm BBMap_38.86.tar.gz

# Add BBTools to the PATH
ENV PATH="/opt/bbmap:${PATH}"
COPY scripts /opt/scripts
RUN chmod a+x /opt/scripts

# Set a default command
CMD ["bbmap.sh"]