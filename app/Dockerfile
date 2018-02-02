# if you're doing anything beyond your local machine, please pin this to a specific version at https://hub.docker.com/_/node/
FROM debian

RUN apt-get update && apt-get install -y wget bzip2

RUN mkdir -p /opt/app

# install dependencies first, in a different location for easier app bind mounting for local development
WORKDIR /opt

# install python stack
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda
ENV PATH /opt/miniconda/bin:$PATH

RUN pip install numpy scipy

RUN apt-get update
RUN apt-get -y install libsuitesparse-dev
RUN pip install cython
RUN pip install -v scikit-sparse
RUN pip install websockets
RUN pip install Flask Flask-Sockets gevent gevent-websocket greenlet gunicorn itsdangerous Jinja2 MarkupSafe redis Werkzeug 
# copy in our source code last, as it changes the most
WORKDIR /opt/app
COPY . /opt/app

CMD [ "python", "-u", "main.py" ]