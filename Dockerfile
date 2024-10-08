FROM fredhutch/r-shiny-server-base:4.4.1
RUN apt-get update -y && apt-get install -y pandoc libudunits2-dev libproj22 libgdal-dev patch
RUN R -q -e 'install.packages(c("shiny", "Seurat", "ggplot2", "plotly", "glue", "reactable", "hdf5r"))'
RUN R -q -e 'remove.packages("bslib")'
RUN R -q -e 'install.packages("bslib", type="source", repos="https://cran.r-project.org")'


RUN rm -rf /srv/shiny-server/

ADD ./app /srv/shiny-server/
ADD check.R /tmp/

RUN chown -R shiny:shiny /srv/shiny-server/

EXPOSE 3838

WORKDIR /srv/shiny-server

RUN R -f /tmp/check.R --args shiny Seurat ggplot2 plotly bslib glue reactable hdf5r

RUN rm /tmp/check.R

ENV SHINY_LOG_STDERR=1 

CMD ["/usr/bin/shiny-server"]
