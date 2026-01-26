CC = g++
CFLAGS = -Wall -I/usr/include/ -I/usr/include/eigen3 -I/usr/include/eigen3 -L/usr/local/lib
LIBS = -lgsl -lgslcblas -lm -larmadillo

BINS = bin/crs bin/delta2tau bin/distribute bin/eV2ex bin/eV2nm bin/eV2rads \
	 bin/Esat bin/fro bin/Gap bin/lyns bin/nom bin/nsISS bin/nsh bin/nsn \
	 bin/nsx bin/nss bin/oap bin/oGp bin/rho2ome_sp bin/sap bin/sfr bin/sgl \
	 bin/sGp bin/two_zero bin/vrb

all: $(BINS)

bin:
	mkdir -p bin

bin/crs: src/cross_section.cxx | bin
	$(CC) $(CFLAGS) src/cross_section.cxx -o bin/crs $(LIBS)

bin/delta2tau: src/delta2tau.cxx | bin
	$(CC) src/delta2tau.cxx -o bin/delta2tau

bin/distribute: src/distribute.cxx | bin
	$(CC) src/distribute.cxx -o bin/distribute

bin/eV2ex: src/eV2ex.cxx | bin
	$(CC) $(CFLAGS) src/eV2ex.cxx -o bin/eV2ex $(LIBS)

bin/eV2nm: src/eV2nm.cxx | bin
	$(CC) $(CFLAGS) src/eV2nm.cxx -o bin/eV2nm $(LIBS)

bin/eV2rads: src/eV2rads.cxx | bin
	$(CC) src/eV2rads.cxx -o bin/eV2rads

bin/Esat: src/Esat.cxx | bin
	$(CC) $(CFLAGS) src/Esat.cxx -o bin/Esat $(LIBS)

bin/fro: src/frohlich.cxx | bin
	$(CC) $(CFLAGS) src/frohlich.cxx -o bin/fro $(LIBS)

bin/Gap: src/nanoshell_G_p3.cxx | bin
	$(CC) $(CFLAGS) src/nanoshell_G_p3.cxx -o bin/Gap $(LIBS)

bin/lyns: src/lycurguseV_ns.cxx | bin
	$(CC) $(CFLAGS) src/lycurguseV_ns.cxx -o bin/lyns $(LIBS)

bin/nom: src/numOme_gsl.cxx | bin
	$(CC) $(CFLAGS) src/numOme_gsl.cxx -o bin/nom $(LIBS)

bin/nsISS: src/nanoshell_ISS.cxx | bin
	$(CC) $(CFLAGS) src/nanoshell_ISS.cxx -o bin/nsISS $(LIBS)

bin/nsh: src/nanoshell.cxx | bin
	$(CC) $(CFLAGS) src/nanoshell.cxx -o bin/nsh $(LIBS)

bin/nsn: src/nanoshell_num.cxx | bin
	$(CC) $(CFLAGS) src/nanoshell_num.cxx -o bin/nsn $(LIBS)

bin/nsx: src/nanoshell_num_anl.cxx | bin
	$(CC) $(CFLAGS) src/nanoshell_num_anl.cxx -o bin/nsx $(LIBS)

bin/nss: src/nanoshell_ss_spe.cxx | bin
	$(CC) $(CFLAGS) src/nanoshell_ss_spe.cxx -o bin/nss $(LIBS)

bin/oap: src/nanoshell_ome_al_p3.cxx | bin
	$(CC) $(CFLAGS) src/nanoshell_ome_al_p3.cxx -o bin/oap $(LIBS)

bin/oGp: src/nanoshell_omeG_p3.cxx | bin
	$(CC) $(CFLAGS) src/nanoshell_omeG_p3.cxx -o bin/oGp $(LIBS)

bin/rho2ome_sp: src/rho2ome_sp.cxx | bin
	$(CC) $(CFLAGS) src/rho2ome_sp.cxx -o bin/rho2ome_sp $(LIBS)

bin/sap: src/single_ome_al_p3.cxx | bin
	$(CC) $(CFLAGS) src/single_ome_al_p3.cxx -o bin/sap $(LIBS)

bin/sfr: src/sfrohlich.cxx | bin
	$(CC) $(CFLAGS) src/sfrohlich.cxx -o bin/sfr $(LIBS)

bin/sgl: src/single.cxx | bin
	$(CC) $(CFLAGS) src/single.cxx -o bin/sgl $(LIBS)

bin/sGp: src/single_omeG_p3.cxx | bin
	$(CC) $(CFLAGS) src/single_omeG_p3.cxx -o bin/sGp $(LIBS)

bin/two_zero: src/two_zero.cxx | bin
	$(CC) src/two_zero.cxx -o bin/two_zero

bin/vrb: src/variables.cxx | bin
	$(CC) $(CFLAGS) src/variables.cxx -o bin/vrb $(LIBS)

clean:
	rm -f $(BINS)