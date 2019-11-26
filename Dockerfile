# tryFBA1
FROM publicisworldwide/python-conda
MAINTAINER  yangyi@tib.cas.cn
COPY tryFBA1/ /home/tryFBA1
WORKDIR /home/tryFBA1/
RUN pip install -r requirements.txt
CMD python tryFBA1.py
