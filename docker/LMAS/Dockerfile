FROM python:3.8.0-slim-buster
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

WORKDIR /data

RUN apt-get update && apt-get install -y build-essential procps
RUN pip install --upgrade pip setuptools wheel
#RUN apt-get install libxml2-dev libxmlsec1-dev

# Install requirements
COPY requirements.txt ./
RUN pip install  --no-cache-dir --use-feature=2020-resolver -r requirements.txt
