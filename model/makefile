#
#	$Id: Makefile,v 1.1.2.4 1999/06/21 11:42:12 severijn Exp $
#

INCS =	common_blocks.h dimensions.h stmt_funcs.h

OBJS =	micom1.o  barotp.o blkdat.o cnuity.o cycl2.o \
	dpthuv.o dpudpv.o geopar.o restin.o restout.o \
	momeq1.o momeq2.o avgout.o sshout.o \
	tsadvc.o boundp.o boundu.o bounduu.o \
        boundv.o boundpb.o boundub.o boundvb.o boundvv.o \
        extra_output.o cyclo.o vort1.o sshini.o \
        boundpne1.o

SCRIPTS=compare_script

SRCS =	micom_ssd.f barotp.f blkdat.f cnuity.f cycl2.f \
	dpthuv.f dpudpv.f geopar.f restin.f restout.f \
	momeq1.f momeq2.f avgout.f sshout.f \
	tsadvc.f boundp.f boundu.f bounduu.f \
        boundv.f boundpb.f boundub.f boundvb.f boundvv.f \
        extra_output.f cyclo.f vort1.f sshini.f \
        boundpne1.f
BUILDDIR=../run
CMD=$(BUILDDIR)/micom1.x
FC=gfortran
LD=$(FC)
# Portland Group pgf90 flags
# FFLAGS=-g -C -Mpreprocess -i4 -r8 -Ktrap=fp
# GCC gfortran flags
FFLAGS= -g -O2 -fbounds-check -std=legacy -fdefault-real-8 -fdefault-double-8 -ffpe-trap=invalid,zero,overflow

all: $(CMD)
clean:
	rm -f $(MICOM) $(OBJS)

demo: $(CMD)

tar: $(MAKEFILE) $(SCRIPTS) $(INCS) $(SRCS)
	@ls $(MAKEFILE) $(SCRIPTS) $(INCS) $(SRCS) > $(TARFILES)

$(CMD): $(OBJS)
	$(LD) $(LDFLAGS) -o $(CMD) $(OBJS)

micom1.o: micom_ssd.f dimensions.h common_blocks.h stmt_funcs.h
	$(FC) $(FFLAGS) -c micom_ssd.f -o micom1.o
barotp.o: barotp.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c barotp.f
blkdat.o: blkdat.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c blkdat.f
cnuity.o: cnuity.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c cnuity.f
dpthuv.o: dpthuv.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c dpthuv.f
dpudpv.o: dpudpv.f dimensions.h
	$(FC) $(FFLAGS) -c dpudpv.f
geopar.o: geopar.f dimensions.h common_blocks.h stmt_funcs.h
	$(FC) $(FFLAGS) -c geopar.f
momeq1.o: momeq1.f dimensions.h common_blocks.h stmt_funcs.h
	$(FC) $(FFLAGS) -c momeq1.f
momeq2.o: momeq2.f dimensions.h common_blocks.h stmt_funcs.h
	$(FC) $(FFLAGS) -c momeq2.f
tsadvc.o: tsadvc.f dimensions.h common_blocks.h stmt_funcs.h
	$(FC) $(FFLAGS) -c tsadvc.f
cyclo.o: cyclo.f dimensions.h common_blocks.h stmt_funcs.h
	$(FC) $(FFLAGS) -c cyclo.f
cycl2.o: cycl2.f dimensions.h common_blocks.h stmt_funcs.h
	$(FC) $(FFLAGS) -c cycl2.f
boundp.o: boundp.f dimensions.h
	$(FC) $(FFLAGS) -c boundp.f
boundu.o: boundu.f dimensions.h
	$(FC) $(FFLAGS) -c boundu.f
boundv.o: boundv.f dimensions.h
	$(FC) $(FFLAGS) -c boundv.f
bounduu.o: bounduu.f dimensions.h
	$(FC) $(FFLAGS) -c bounduu.f
boundvv.o: boundvv.f dimensions.h
	$(FC) $(FFLAGS) -c boundvv.f
boundpb.o: boundpb.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c boundpb.f
boundpne1.o: boundpne1.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c boundpne1.f
boundub.o: boundub.f dimensions.h
	$(FC) $(FFLAGS) -c boundub.f
boundvb.o: boundvb.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c boundvb.f
vort1.o: vort1.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c vort1.f
extra_output.o: extra_output.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c extra_output.f
avgout.o: avgout.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c avgout.f
sshout.o: sshout.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c sshout.f
sshini.o: sshini.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c sshini.f
restin.o: restin.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c restin.f
restout.o: restout.f dimensions.h common_blocks.h
	$(FC) $(FFLAGS) -c restout.f
