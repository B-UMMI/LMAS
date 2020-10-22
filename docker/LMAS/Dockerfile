FROM python:3.9
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

WORKDIR /data

RUN pip install --upgrade pip setuptools wheel

# Install requirements
COPY requirements.txt ./
RUN pip install  --no-cache-dir -r requirements.txt