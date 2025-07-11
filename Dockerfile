FROM python:3.10.16-bullseye

RUN mkdir /workspace
COPY requirements.txt /workspace
WORKDIR /workspace

RUN apt-get update && apt-get install -y \
    sudo \
    tree \
    man  

RUN pip install -r requirements.txt

RUN groupadd -g 1012 admin
RUN useradd admin -u 1012 -g 1012 -m -s /bin/bash && echo "admin:admin" | chpasswd && adduser admin sudo
USER admin

# Install the image
# docker build -f Dockerfile -t uslenghi_image .

# Run the container with the following command:
# docker run -it --name uslenghi_container -v /mnt/ABLAZIONE/raw_data:/workspace/raw_data:ro -v /dumpall/uslenghi/LA_entropy/processed_data:/workspace/processed_data -v /home/uslenghi/LA_entropy:/workspace uslenghi_image bash &
