# syntax= docker/dockerfile:1

FROM python:3.9.9

LABEL author="Priya Kempanna"
LABEL author="Lara Schmalenstroer"
LABEL author="Rohitha Ravinder"
LABEL author="Shubhi Ambast"
LABEL author_email="kpriyaa97@gmail.com"
LABEL author_email="lara.schmalenstroer@gmail.com"
LABEL author_email="rohitha0112@gmail.com"
LABEL author_email="shubhiambast@gmail.com"
LABEL description="Group03 Protein Identification Mass Spectrometry Package"


WORKDIR Desktop/group_3/

COPY . .

RUN pip3 install -e mass_spectrum

ENV FLASK_PORT=5000
EXPOSE ${FLASK_PORT}

ENTRYPOINT ["python3", "frontend/run.py"]
