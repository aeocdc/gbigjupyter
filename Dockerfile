# gBIGjupyter enviroments
FROM python
MAINTAINER  yangyi@tib.cas.cn
#USER root
RUN apt-get update && apt-get install -y ocl-icd-opencl-dev
RUN wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3
RUN /opt/miniconda3/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ && \
    /opt/miniconda3/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/ && \
    /opt/miniconda3/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/mro/ && \
    /opt/miniconda3/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
    /opt/miniconda3/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
    /opt/miniconda3/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
RUN /opt/miniconda3/bin/conda install -y cas-offinder
RUN pip install jupyter
RUN pip install pocl
#COPY requirement.txt ./requirement.txt
#RUN pip install esdk-obs-python --trusted-host pypi.org
#RUN pip install -r requirement.txt
#ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN mkdir /home/gbig/
WORKDIR /home/gbig
COPY ref ./ref
COPY computational_tools ./computational_tools
COPY python_script ./python_script
COPY input ./input
CMD jupyter-notebook --ip=0.0.0.0 --no-browser --allow-root
