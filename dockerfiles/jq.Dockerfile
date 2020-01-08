FROM stedolan/jq

RUN apt-get update && apt-get install -y gawk \
    && apt-get autoremove -y