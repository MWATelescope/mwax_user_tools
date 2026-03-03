FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PSRHOME=/home/psr/software
ENV OSTYPE=linux

# Create psr user
RUN useradd -m -s /bin/bash psr && mkdir -p $PSRHOME && chown -R psr:psr $PSRHOME

# ── System dependencies ────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y \
    build-essential autoconf autotools-dev automake autogen libtool \
    pkg-config cmake csh \
    gcc gfortran g++ \
    wget git curl \
    libcfitsio-dev \
    libfftw3-dev libfftw3-bin libfftw3-single3 \
    libpng-dev libpnglite-dev \
    libx11-dev \
    libhdf5-dev libhdf5-serial-dev \
    libxml2-dev \
    libltdl-dev \
    libgsl-dev gsl-bin \
    pgplot5 \
    python3 python3-dev python3-pip \
    swig \
    && rm -rf /var/lib/apt/lists/*

ENV PGPLOT_DIR=/usr/lib/pgplot5
ENV PGPLOT_FONT=/usr/lib/pgplot5/grfont.dat
ENV PGPLOT_DEV=/xs

USER psr
WORKDIR $PSRHOME

# ── 1. calceph ────────────────────────────────────────────────────────────────
RUN wget -q https://www.imcce.fr/content/medias/recherche/equipes/asd/calceph/calceph-3.5.1.tar.gz && \
    tar -xf calceph-3.5.1.tar.gz && mv calceph-3.5.1 calceph && rm calceph-3.5.1.tar.gz
WORKDIR $PSRHOME/calceph
RUN ./configure --prefix=$PSRHOME/calceph/install \
                --enable-shared --enable-static --with-pic \
                --enable-fortran --enable-thread && \
    make -j$(nproc) && make install

# ── 2. TEMPO2 ─────────────────────────────────────────────────────────────────
RUN git clone https://bitbucket.org/psrsoft/tempo2.git $PSRHOME/tempo2-src
WORKDIR $PSRHOME/tempo2-src
RUN sync && perl -pi -e 's/chmod \+x/#chmod +x/' bootstrap && \
    ./bootstrap && \
    ./configure --prefix=$PSRHOME/tempo2/T2runtime \
                --enable-shared --enable-static --with-pic \
                --with-calceph=$PSRHOME/calceph/install \
                CPPFLAGS="-I$PSRHOME/calceph/install/include" \
                LDFLAGS="-L$PSRHOME/calceph/install/lib" \
                LIBS="-lcalceph" \
                F77=gfortran \
                FFLAGS="-fallow-argument-mismatch" && \
    make -j$(nproc) && make install && \
    make plugins && \
    make plugins-install DESTDIR=$PSRHOME/tempo2/T2runtime plugindir=$PSRHOME/tempo2/T2runtime/plugins || \
    echo "WARNING: plugins-install failed (non-fatal, skipping)" && \
    cp *.h $PSRHOME/tempo2/T2runtime/include/ 2>/dev/null || true

# ── 3. psrcat ─────────────────────────────────────────────────────────────────
RUN wget -q https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz && \
    tar -xf psrcat_pkg.tar.gz -C $PSRHOME && rm psrcat_pkg.tar.gz
WORKDIR $PSRHOME/psrcat_tar
RUN /bin/bash makeit

# ── 4. PSRCHIVE ───────────────────────────────────────────────────────────────
RUN git clone git://git.code.sf.net/p/psrchive/code $PSRHOME/psrchive
WORKDIR $PSRHOME/psrchive
RUN ./bootstrap && \
    ./configure --prefix=$PSRHOME/psrchive/install \
                --enable-shared --enable-static \
                F77=gfortran && \
    make -j$(nproc) && make install
WORKDIR /home/psr
RUN echo "Predictor::default = tempo2" >> .psrchive.cfg && \
    echo "Predictor::policy = default" >> .psrchive.cfg

# ── 5. DSPSR with VDIF ────────────────────────────────────────────────────────
# All paths set together in one ENV block to avoid the cascading
# $VAR=$VAR:... expansion problem that silently empties variables in Docker.
# PKG_CONFIG_PATH is critical — dspsr configure uses pkg-config to find psrchive.
ENV PATH=$PSRHOME/calceph/install/bin:\
$PSRHOME/tempo2/T2runtime/bin:\
$PSRHOME/psrcat_tar:\
$PSRHOME/psrchive/install/bin:\
/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    LD_LIBRARY_PATH=$PSRHOME/calceph/install/lib:\
$PSRHOME/tempo2/T2runtime/lib:\
$PSRHOME/psrchive/install/lib:\
$PSRHOME/psrdada/install/lib \
    C_INCLUDE_PATH=$PSRHOME/calceph/install/include:\
$PSRHOME/tempo2/T2runtime/include:\
$PSRHOME/psrchive/install/include \
    PKG_CONFIG_PATH=$PSRHOME/psrchive/install/lib/pkgconfig:\
$PSRHOME/calceph/install/lib/pkgconfig \
    TEMPO2=$PSRHOME/tempo2/T2runtime \
    PSRCAT_FILE=$PSRHOME/psrcat_tar/psrcat.db \
    PSRDADA_HOME=$PSRHOME/psrdada/install

# ── 5. PSRDADA ────────────────────────────────────────────────────────────────
RUN git clone git://git.code.sf.net/p/psrdada/code $PSRHOME/psrdada
WORKDIR $PSRHOME/psrdada
RUN ./bootstrap && \
    ./configure --prefix=$PSRHOME/psrdada/install && \
    make -j$(nproc) && make install

# ── 6. DSPSR with VDIF + DADA ─────────────────────────────────────────────────
RUN git clone git://git.code.sf.net/p/dspsr/code $PSRHOME/dspsr
WORKDIR $PSRHOME/dspsr
RUN echo "Checking for Tempo2::Generator in psrchive libs:" && \
    nm -D $PSRHOME/psrchive/install/lib/libpsrmore.so 2>/dev/null | grep -i "T2Generator\|Tempo2.*Generator" || \
    nm -D $PSRHOME/psrchive/install/lib/libpsrutil.so 2>/dev/null | grep -i "T2Generator\|Tempo2.*Generator" || \
    nm -D $PSRHOME/psrchive/install/lib/libpsrbase.so 2>/dev/null | grep -i "T2Generator\|Tempo2.*Generator" || \
    echo "NOT FOUND in any psrchive lib — psrchive needs --with-tempo2 at configure time"

# Rebuild psrchive with explicit tempo2 support
WORKDIR $PSRHOME/psrchive
RUN ./configure --prefix=$PSRHOME/psrchive/install \
                --enable-shared --enable-static \
                --with-tempo2-dir=$PSRHOME/tempo2/T2runtime \
                F77=gfortran \
                CPPFLAGS="-I$PSRHOME/tempo2/T2runtime/include \
                           -I$PSRHOME/tempo2-src" \
                LDFLAGS="-L$PSRHOME/tempo2/T2runtime/lib" && \
    make -j$(nproc) && make install

WORKDIR $PSRHOME/dspsr
RUN echo "vdif dada fits sigproc" > backends.list && \
    ./bootstrap && \
    ./configure --prefix=$PSRHOME/dspsr/install \
                F77=gfortran \
                CPPFLAGS="-I$PSRHOME/psrchive/install/include \
                           -I$PSRHOME/psrchive/install/include/epsic \
                           -I$PSRHOME/psrchive/Util/tempo2 \
                           -I$PSRHOME/calceph/install/include \
                           -I$PSRHOME/tempo2/T2runtime/include \
                           -I$PSRHOME/tempo2-src \
                           -I$PSRHOME/psrdada/install/include" \
                LDFLAGS="-L$PSRHOME/psrchive/install/lib \
                          -L$PSRHOME/calceph/install/lib \
                          -L$PSRHOME/tempo2/T2runtime/lib \
                          -L$PSRHOME/psrdada/install/lib" \
                LIBS="-lpsrdada -lpsrmore -lpsrbase -lpsrutil" && \
    make -j$(nproc) && make install

ENV PATH=$PATH:$PSRHOME/dspsr/install/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PSRHOME/dspsr/install/lib

WORKDIR /home/psr