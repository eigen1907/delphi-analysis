# Build with --build-context generator=/path/to/delphisimulation/pythia8.
# DELSIM/DELANA come unchanged from this public x86 image.
FROM docker.io/jingyucms/delphi-pythia8:v2.4@sha256:501ae4d22a2a4bbcc76f9d8d12536377b8e06c109caa35da6081310aed14ecc3

USER root
RUN dnf -y install epel-release dnf-plugins-core \
    && dnf -y config-manager --set-enabled crb \
    && dnf -y install cmake gcc-gfortran make git \
        root root-core root-io root-tree root-hist root-graf root-graf3d \
        root-mathcore root-matrix root-physics root-genvector root-smatrix \
        root-montecarlo-eg HepMC HepMC-devel HepMC3 HepMC3-devel \
    && dnf clean all

# The generator uses PHOTOS++ through HepMC2; build both supported interfaces.
RUN git init /tmp/photospp \
    && git -C /tmp/photospp fetch --depth 1 https://gitlab.cern.ch/photospp/photospp.git \
        9f860e003c9aa9898e02dc3fcf78ffb2ca81aa17 \
    && git -C /tmp/photospp checkout --detach FETCH_HEAD \
    && cmake -S /tmp/photospp -B /tmp/photospp/build \
        -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release \
        -DPhotospp_ENABLE_HEPMC2=ON -DPhotospp_ENABLE_HEPMC3=ON \
        -DPhotospp_ENABLE_EXAMPLES=OFF -DPhotospp_ENABLE_MC-TESTER=OFF \
        -DPhotospp_ENABLE_DOXYGEN=OFF -DPhotospp_ENABLE_LATEX=OFF \
        -DPYTHIA8_DIR=/usr -DHepMC3_DIR=/usr \
    && cmake --build /tmp/photospp/build --parallel 2 \
    && cmake --install /tmp/photospp/build \
    && rm -rf /tmp/photospp

# The native NanoAOD converter links FastJet, which is absent from EPEL 9.
RUN curl -fL https://fastjet.fr/repo/fastjet-3.4.3.tar.gz -o /tmp/fastjet.tar.gz \
    && tar -xzf /tmp/fastjet.tar.gz -C /tmp \
    && cd /tmp/fastjet-3.4.3 \
    && ./configure --prefix=/usr/local --enable-shared --enable-allcxxplugins \
    && make -j2 && make install \
    && rm -rf /tmp/fastjet-3.4.3 /tmp/fastjet.tar.gz

ENV Photospp_DIR=/usr/local
RUN printf '/usr/local/lib\n/usr/local/lib64\n' > /etc/ld.so.conf.d/photospp.conf \
    && ldconfig

# Import the upstream generator read-only at build time; keep binaries outside
# /work so mounting an event-production directory cannot hide them.
COPY --from=generator Makefile pythia8_generate.cpp /opt/generator/
# Build/link helpers and PHDST runtime locations.
ENV GROUP_DIR="/delphi/releases/almalinux-9-x86_64/latest"
ENV CERN="${GROUP_DIR}/cern" \
    CERN_ROOT="${GROUP_DIR}/cern/pro" \
    CERN_LIB="${GROUP_DIR}/cern/pro/lib" \
    CERN_PAM="${GROUP_DIR}/cern/pro/src/car" \
    DELPHI="${GROUP_DIR}/dstana/161018" \
    DELPHI_LIB="${GROUP_DIR}/dstana/161018/lib" \
    DELPHI_PAM="${GROUP_DIR}/dstana/161018/src/car" \
    DELPHI_CRA="${GROUP_DIR}/dstana/161018/src/car" \
    DELPHI_BLKD="${GROUP_DIR}/dstana/161018/blkd" \
    DELPHI_DAT="${GROUP_DIR}/dstana/161018/dat" \
    DELPHI_BIN="${GROUP_DIR}/bin" \
    ADDLIB="${GROUP_DIR}/dstana/161018/blkd/delblkd.o" \
    PLINAM="UNIX" \
    DELSIM_ROOT="${GROUP_DIR}/simana" \
    DELPHI_RELEASE="latest" \
    DELPHI_PATH="${GROUP_DIR}/bin" \
    GROUPPATH="${GROUP_DIR}/bin" \
    DELPHI_DDB_BIN="${GROUP_DIR}/ddb/bin" \
    DELPHI_DDB_DIR="${GROUP_DIR}/ddb" \
    DELPHI_INSTALL_DIR="/delphi" \
    DELPHI_DATA_ROOT="/eos/opendata/delphi" \
    DELPHI_DDB="/eos/opendata/delphi/condition-data" \
    DELPHI_ZIP="on" \
    DELPHI_XRD="" \
    DELPHI_BATCAVE="/delphi/pdl" \
    DES_HOME="${GROUP_DIR}/evserv" \
    GRA_PLACE="${GROUP_DIR}/delgra" \
    OPGS_RUN_TIME="${GROUP_DIR}/delgra/run_time" \
    OPENPHIGS="${GROUP_DIR}/openphigs" \
    PATH="${GROUP_DIR}/scripts:/delphi/scripts:${GROUP_DIR}/bin:${GROUP_DIR}/cern/pro/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    PYTHIA8DATA="/usr/share/Pythia8/xmldoc" \
    LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/usr/lib64:/lib64"
RUN make -C /opt/generator pythia8_generate \
    && cp /opt/generator/pythia8_generate /usr/local/bin/pythia8_generate

WORKDIR /work
ENTRYPOINT ["/bin/bash"]
