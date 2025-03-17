FROM python:3.10.16-bullseye

RUN mkdir /workspace
COPY requirements.txt /workspace
WORKDIR /workspace

RUN apt-get update && apt-get install -y \
    tree

RUN pip install -r requirements.txt

# Install the image
# docker build -f Dockerfile -t uslenghi_image .

# Run the container with the following command:
# docker run -it --name uslenghi_container -v /mnt/ABLAZIONE/raw_data:/workspace/raw_data:ro  -v /home/uslenghi/LA_entropy:/workspace uslenghi_image bash &
