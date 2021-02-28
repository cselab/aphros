$(WRK)/.dir:
	mkdir -p $(WRK)/color && \
	mkdir -p $(WRK)/distr && \
	mkdir -p $(WRK)/dump && \
	mkdir -p $(WRK)/func && \
	mkdir -p $(WRK)/geom && \
	mkdir -p $(WRK)/inside && \
	mkdir -p $(WRK)/linear && \
	mkdir -p $(WRK)/march && \
	mkdir -p $(WRK)/parse && \
	mkdir -p $(WRK)/solver && \
	mkdir -p $(WRK)/util && \
	mkdir -p $(WRK)/young && \
	touch $@
